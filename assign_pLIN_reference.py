#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Assign pLIN codes to ALL reference plasmid sequences (~72,959).

Pipeline:
  Phase 1 — Compute 4-mer frequency vectors for all reference sequences
  Phase 2 — Classify Inc types using the pre-trained KNN classifier
  Phase 3 — Per-Inc-group single-linkage clustering → pLIN codes
  Phase 4 — Output combined table

Usage:
  python assign_pLIN_reference.py [--resume] [--max-group-size 25000]
"""

import os
import sys
import glob
import json
import time
import hashlib
import argparse
import subprocess
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from itertools import product as iter_product
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.neighbors import KNeighborsClassifier
from Bio import SeqIO

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CLASSIFIER_PATH = os.path.join(BASE_DIR, "data", "inc_classifier.npz")
SEQUENCES_FASTA = os.path.join(BASE_DIR, "sequences.fasta")
REFERENCE_DIR = os.path.join(BASE_DIR, "reference")
TRAINING_DIR = os.path.join(BASE_DIR, "plasmid_sequences_for_training")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")

# Checkpoints
VECTORS_CHECKPOINT = os.path.join(OUTPUT_DIR, "reference_kmer_vectors.npz")
INC_CHECKPOINT = os.path.join(OUTPUT_DIR, "reference_inc_classifications.tsv")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "pLIN_reference_assignments.tsv")

# pLIN hierarchical thresholds
PLIN_THRESHOLDS = {
    "A": 0.150,   # ~85% ANI — L1 (Family)
    "B": 0.100,   # ~90% ANI — L2 (Subfamily)
    "C": 0.050,   # ~95% ANI — L3 (Cluster)
    "D": 0.020,   # ~98% ANI — L4 (Subcluster)
    "E": 0.010,   # ~99% ANI — L5 (Clone)
    "F": 0.001,   # ~99.9% ANI — L6 (Strain)
}

# Inc classification thresholds (from plin_app.py)
INC_CONFIDENCE_THRESHOLD = 0.40
MULTI_INC_THRESHOLD = 0.15

# See assign_pLIN.py for rationale: legacy training-provenance labels are
# mapped onto current PlasmidFinder nomenclature for reporting.
RESOLVED_INC_TYPE = {
    "IncAC2": "IncA/IncC",
}


def resolve_inc_type(inc_type):
    return RESOLVED_INC_TYPE.get(inc_type, inc_type)


def run_version_stamp(training_fingerprint=None):
    """Reproducibility stamp — see assign_pLIN.py::run_version_stamp for rationale."""
    try:
        git_hash = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=BASE_DIR,
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        git_hash = "unknown"

    return {
        "run_timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "git_commit": git_hash,
        "training_input_fingerprint": training_fingerprint,
        "classifier_path": CLASSIFIER_PATH,
    }

# 4-mer setup
K = 4
BASES = "ACGT"
KMERS = ["".join(p) for p in iter_product(BASES, repeat=K)]
KMER_INDEX = {km: i for i, km in enumerate(KMERS)}
N_KMERS = len(KMERS)  # 256


def stable_cluster_ids(raw_labels, member_keys):
    """Renumber scipy fcluster labels into a deterministic scheme.

    scipy.cluster.hierarchy.fcluster assigns cluster-ID integers based on
    internal linkage-tree traversal order, which is not guaranteed stable
    across runs even when cluster MEMBERSHIP is identical. Renumbering by
    each cluster's lexicographically smallest member key (plasmid accession)
    makes cluster IDs reproducible across runs and pipeline versions.

    NOTE: this alone only guarantees reproducibility for a FIXED input set.
    Adding new sequences to the database can change which member is
    lexicographically smallest within a cluster, which reassigns the ID even
    when the underlying cluster relationships are unchanged. For stability
    across database growth, see `frozen_seed_cluster_ids` below, which is
    the scheme actually used once a frozen training-set numbering exists.
    """
    raw_labels = np.asarray(raw_labels)
    cluster_min_key = {}
    for raw_id, key in zip(raw_labels, member_keys):
        if raw_id not in cluster_min_key or key < cluster_min_key[raw_id]:
            cluster_min_key[raw_id] = key

    ordered_raw_ids = sorted(cluster_min_key, key=lambda rid: cluster_min_key[rid])
    remap = {raw_id: new_id for new_id, raw_id in enumerate(ordered_raw_ids, start=1)}
    return np.array([remap[rid] for rid in raw_labels])


def frozen_seed_cluster_ids(raw_labels, member_keys, frozen_id_of_key):
    """Renumber clusters so that IDs inherited from a frozen prior run are
    preserved, and only genuinely new clusters receive fresh IDs.

    `frozen_id_of_key`: dict mapping a subset of `member_keys` (the ones
    that were already assigned a canonical ID in a prior, frozen run — e.g.
    the training-only pLIN_assignments.tsv) to that canonical integer ID.

    For each cluster found in THIS run, if any of its members carry a
    frozen ID, the cluster inherits that ID (majority vote if it somehow
    spans more than one frozen ID — this should be rare and only happens
    when growth merges two previously-separate frozen clusters, which we
    resolve by keeping the smaller original frozen ID for continuity).
    Clusters with no frozen members at all (i.e. built entirely from newly
    added reference sequences) get fresh sequential IDs continuing after
    the maximum frozen ID, keyed by lexicographically smallest member for
    reproducibility among themselves.
    """
    raw_labels = np.asarray(raw_labels)

    # Group member indices by raw cluster label
    cluster_members = {}
    for i, (raw_id, key) in enumerate(zip(raw_labels, member_keys)):
        cluster_members.setdefault(raw_id, []).append((i, key))

    max_frozen_id = max(frozen_id_of_key.values()) if frozen_id_of_key else 0

    remap = {}
    new_clusters = []  # (raw_id, min_key) for clusters with no frozen anchor
    for raw_id, members in cluster_members.items():
        frozen_ids_present = sorted(set(
            frozen_id_of_key[key] for _, key in members if key in frozen_id_of_key
        ))
        if frozen_ids_present:
            # Inherit the smallest frozen ID seen (stable, deterministic
            # tie-break if growth merged two previously-distinct clusters)
            remap[raw_id] = frozen_ids_present[0]
        else:
            min_key = min(key for _, key in members)
            new_clusters.append((raw_id, min_key))

    # Assign fresh IDs to genuinely new clusters, ordered by min member key
    # for reproducibility, continuing the numbering after the frozen max.
    new_clusters.sort(key=lambda t: t[1])
    for offset, (raw_id, _) in enumerate(new_clusters, start=1):
        remap[raw_id] = max_frozen_id + offset

    return np.array([remap[rid] for rid in raw_labels])


# ── Phase 1: Load sequences & compute 4-mer vectors ──────────────────────────

def kmer_vector(sequence):
    """Compute normalised 4-mer frequency vector (sliding window, fast)."""
    seq = sequence.upper()
    counts = np.zeros(N_KMERS, dtype=np.float64)
    for i in range(len(seq) - K + 1):
        kmer = seq[i:i + K]
        if kmer in KMER_INDEX:
            counts[KMER_INDEX[kmer]] += 1
    total = counts.sum()
    if total > 0:
        counts /= total
    return counts


def load_reference_sequences():
    """Load all reference sequences from sequences.fasta or reference/ dir."""
    records = []

    if os.path.exists(SEQUENCES_FASTA):
        print(f"  Reading from {SEQUENCES_FASTA} ...", flush=True)
        for rec in SeqIO.parse(SEQUENCES_FASTA, "fasta"):
            records.append({
                "plasmid_id": rec.id,
                "sequence": str(rec.seq),
                "length": len(rec.seq),
                "source": "reference",
            })
            if len(records) % 10000 == 0:
                print(f"    {len(records):,} sequences loaded ...", flush=True)
    elif os.path.isdir(REFERENCE_DIR):
        print(f"  Reading from {REFERENCE_DIR}/ ...", flush=True)
        fasta_files = sorted(glob.glob(os.path.join(REFERENCE_DIR, "*.fasta")))
        for fpath in fasta_files:
            for rec in SeqIO.parse(fpath, "fasta"):
                records.append({
                    "plasmid_id": rec.id,
                    "sequence": str(rec.seq),
                    "length": len(rec.seq),
                    "source": "reference",
                })
            if len(records) % 10000 == 0:
                print(f"    {len(records):,} sequences loaded ...", flush=True)
    else:
        print("ERROR: Neither sequences.fasta nor reference/ directory found.", file=sys.stderr)
        sys.exit(1)

    print(f"  Total reference sequences: {len(records):,}", flush=True)
    return records


def load_training_sequences():
    """Load training sequences (with known Inc types) from training folders."""
    records = []
    for inc_dir in sorted(glob.glob(os.path.join(TRAINING_DIR, "*", "fastas"))):
        inc_name = os.path.basename(os.path.dirname(inc_dir))
        fasta_files = sorted(glob.glob(os.path.join(inc_dir, "*.fasta")))
        for fpath in fasta_files:
            for rec in SeqIO.parse(fpath, "fasta"):
                records.append({
                    "plasmid_id": rec.id,
                    "sequence": str(rec.seq),
                    "length": len(rec.seq),
                    "source": "training",
                    "inc_type": inc_name,
                    "inc_confidence": 1.0,
                    "is_novel": False,
                    "is_multiple": False,
                    "inc_secondary": "",
                })
    print(f"  Training sequences: {len(records):,} across "
          f"{len(set(r['inc_type'] for r in records))} Inc groups", flush=True)
    return records


def compute_vectors(records, checkpoint_path, resume=False):
    """Compute 4-mer vectors for all records, with checkpoint support."""
    if resume and os.path.exists(checkpoint_path):
        print(f"  Resuming: loading vectors from {checkpoint_path} ...", flush=True)
        data = np.load(checkpoint_path, allow_pickle=True)
        vectors = data["vectors"]
        ids = list(data["ids"])
        lengths = list(data["lengths"])
        print(f"  Loaded {vectors.shape[0]:,} vectors from checkpoint", flush=True)
        return vectors, ids, lengths

    n = len(records)
    vectors = np.zeros((n, N_KMERS), dtype=np.float32)
    ids = []
    lengths = []

    t0 = time.time()
    for idx, rec in enumerate(records):
        vectors[idx] = kmer_vector(rec["sequence"]).astype(np.float32)
        ids.append(rec["plasmid_id"])
        lengths.append(rec["length"])
        if (idx + 1) % 5000 == 0:
            elapsed = time.time() - t0
            rate = (idx + 1) / elapsed
            eta = (n - idx - 1) / rate
            print(f"    {idx+1:,}/{n:,} vectors computed "
                  f"({rate:.0f}/sec, ETA {eta/60:.1f} min)", flush=True)

    elapsed = time.time() - t0
    print(f"  All {n:,} vectors computed in {elapsed/60:.1f} min", flush=True)

    # Save checkpoint
    np.savez_compressed(checkpoint_path,
                        vectors=vectors,
                        ids=np.array(ids, dtype=object),
                        lengths=np.array(lengths, dtype=np.int64))
    print(f"  Checkpoint saved: {checkpoint_path} "
          f"({os.path.getsize(checkpoint_path)/1024/1024:.1f} MB)", flush=True)

    return vectors, ids, lengths


# ── Phase 2: Classify Inc types ──────────────────────────────────────────────

def classify_inc_types(vectors, ids, lengths, checkpoint_path, resume=False):
    """Classify all sequences by Inc type using KNN."""
    if resume and os.path.exists(checkpoint_path):
        print(f"  Resuming: loading classifications from {checkpoint_path} ...", flush=True)
        df = pd.read_csv(checkpoint_path, sep="\t")
        print(f"  Loaded {len(df):,} classifications from checkpoint", flush=True)
        return df

    print("  Loading KNN classifier ...", flush=True)
    data = np.load(CLASSIFIER_PATH, allow_pickle=True)
    X_train = data["X"]
    y_train = data["y"]
    group_names = [str(g) for g in data["group_names"]]

    knn = KNeighborsClassifier(n_neighbors=5, metric="cosine", weights="distance")
    knn.fit(X_train, y_train)
    print(f"  KNN trained: {len(X_train)} samples, {len(group_names)} groups", flush=True)

    n = len(vectors)
    results = []
    batch_size = 5000
    t0 = time.time()

    for start in range(0, n, batch_size):
        end = min(start + batch_size, n)
        batch = vectors[start:end]
        proba = knn.predict_proba(batch)

        for i in range(len(batch)):
            idx = start + i
            pred_idx = np.argmax(proba[i])
            confidence = float(proba[i][pred_idx])
            best_group = group_names[pred_idx]

            # Check for multiple Inc types
            high_conf = [(group_names[j], float(proba[i][j]))
                         for j in range(len(group_names))
                         if float(proba[i][j]) >= MULTI_INC_THRESHOLD]
            high_conf.sort(key=lambda x: x[1], reverse=True)
            is_multiple = len(high_conf) >= 2
            is_novel = confidence < INC_CONFIDENCE_THRESHOLD

            if is_novel:
                inc_type = "Unknown"
            elif is_multiple:
                inc_type = best_group  # Cluster with primary
            else:
                inc_type = best_group

            inc_secondary = ""
            if is_multiple:
                secondary = [f"{inc}({conf:.2f})" for inc, conf in high_conf[1:3]]
                inc_secondary = "; ".join(secondary)

            results.append({
                "plasmid_id": ids[idx],
                "length_bp": lengths[idx],
                "inc_type": inc_type,
                "inc_confidence": round(confidence, 4),
                "is_novel": is_novel,
                "is_multiple": is_multiple,
                "inc_secondary": inc_secondary,
                "best_match": best_group,
            })

        elapsed = time.time() - t0
        print(f"    {end:,}/{n:,} classified ({elapsed:.0f}s)", flush=True)

    df = pd.DataFrame(results)

    # Save checkpoint
    df.to_csv(checkpoint_path, sep="\t", index=False)
    print(f"  Checkpoint saved: {checkpoint_path}", flush=True)

    # Summary
    print(f"\n  Inc Type Classification Summary:", flush=True)
    print(f"    Total sequences:  {len(df):,}", flush=True)
    print(f"    Classified:       {(~df['is_novel']).sum():,} "
          f"({(~df['is_novel']).mean()*100:.1f}%)", flush=True)
    print(f"    Unknown/Novel:    {df['is_novel'].sum():,} "
          f"({df['is_novel'].mean()*100:.1f}%)", flush=True)
    print(f"    Multiple Inc:     {df['is_multiple'].sum():,}", flush=True)
    print(f"\n  Inc type distribution (classified only):", flush=True)
    classified = df[~df["is_novel"]]
    for inc, count in classified["inc_type"].value_counts().items():
        print(f"    {inc:<12s} {count:>6,}", flush=True)

    return df


# ── Phase 3: Per-Inc-group clustering ────────────────────────────────────────

TRAINING_ASSIGNMENTS_PATH = os.path.join(OUTPUT_DIR, "pLIN_assignments.tsv")


def load_frozen_training_ids(inc_type):
    """Load the frozen (training-only) cluster IDs for one Inc group, per level.

    Returns dict: {bin_name: {plasmid_id: frozen_int_id}}, or {} if the
    frozen training-set assignments file is unavailable (in which case
    cluster_inc_group falls back to the run-local stable_cluster_ids scheme).
    """
    if not hasattr(load_frozen_training_ids, "_cache"):
        if os.path.exists(TRAINING_ASSIGNMENTS_PATH):
            df = pd.read_csv(TRAINING_ASSIGNMENTS_PATH, sep="\t")
            load_frozen_training_ids._cache = df
        else:
            load_frozen_training_ids._cache = None

    df = load_frozen_training_ids._cache
    if df is None:
        return {}

    sub = df[df["inc_type"] == inc_type]
    if sub.empty:
        return {}

    result = {}
    for bname in PLIN_THRESHOLDS:
        col = f"bin_{bname}"
        result[bname] = dict(zip(sub["plasmid_id"], sub[col].astype(int)))
    return result


def cluster_inc_group(vectors, member_keys, max_group_size=25000, inc_type=None):
    """Run single-linkage clustering on vectors within one Inc group.

    When a frozen training-only run exists (output/pLIN_assignments.tsv),
    cluster IDs for this Inc group are anchored to that frozen numbering via
    `frozen_seed_cluster_ids`, so adding new reference sequences does not
    renumber clusters that already existed in the training-only run — only
    genuinely new clusters receive new IDs. Without a frozen reference,
    falls back to `stable_cluster_ids` (deterministic within this run only).
    """
    n = len(vectors)
    if n < 2:
        # Single sequence — assign cluster 1 at all levels
        assignments = {b: np.array([1]) for b in PLIN_THRESHOLDS}
        return assignments

    # Estimate memory
    condensed_size = n * (n - 1) // 2
    mem_gb = condensed_size * 4 / 1e9  # float32

    if n > max_group_size:
        print(f"    WARNING: Group has {n:,} sequences ({mem_gb:.1f} GB). "
              f"Using chunked approach ...", flush=True)

    # Compute pairwise cosine distances
    dist_condensed = pdist(vectors.astype(np.float64), metric="cosine")

    # Single-linkage clustering
    Z = linkage(dist_condensed, method="single")

    frozen = load_frozen_training_ids(inc_type) if inc_type else {}

    # Cut at each threshold, renumbering to a deterministic cluster-ID scheme
    assignments = {}
    for bname, thresh in PLIN_THRESHOLDS.items():
        raw_clusters = fcluster(Z, t=thresh, criterion="distance")
        if bname in frozen and frozen[bname]:
            assignments[bname] = frozen_seed_cluster_ids(raw_clusters, member_keys, frozen[bname])
        else:
            assignments[bname] = stable_cluster_ids(raw_clusters, member_keys)

    return assignments


def phase3_clustering(ref_classifications, ref_vectors, ref_ids,
                      training_records, training_vectors,
                      max_group_size=25000):
    """Per-Inc-group clustering combining training + classified reference sequences."""
    print("\n" + "=" * 70, flush=True)
    print("PHASE 3: Per-Inc-Group Clustering", flush=True)
    print("=" * 70, flush=True)

    # Build lookup from plasmid_id to vector index for reference sequences
    ref_id_to_idx = {pid: i for i, pid in enumerate(ref_ids)}

    # Build lookup for training sequences
    train_id_to_vec = {}
    train_id_to_rec = {}
    for i, rec in enumerate(training_records):
        train_id_to_vec[rec["plasmid_id"]] = training_vectors[i]
        train_id_to_rec[rec["plasmid_id"]] = rec

    # Get all Inc groups (from training + classified reference)
    classified_ref = ref_classifications[~ref_classifications["is_novel"]]
    all_inc_types = sorted(set(
        list(classified_ref["inc_type"].unique()) +
        list(set(r["inc_type"] for r in training_records))
    ))

    # Find overlap: training IDs that also appear in reference
    training_ids = set(r["plasmid_id"] for r in training_records)
    reference_ids = set(ref_ids)
    overlap_ids = training_ids & reference_ids
    if overlap_ids:
        print(f"  Note: {len(overlap_ids):,} sequences appear in both "
              f"training and reference (will use training Inc type)", flush=True)

    all_results = []
    t0 = time.time()

    for inc_type in all_inc_types:
        # Gather training sequences of this Inc type
        train_ids_this = [r["plasmid_id"] for r in training_records
                          if r["inc_type"] == inc_type]
        # Gather classified reference sequences of this Inc type
        ref_rows_this = classified_ref[classified_ref["inc_type"] == inc_type]
        # Exclude reference sequences that are also in training (avoid duplicates)
        ref_ids_this = [pid for pid in ref_rows_this["plasmid_id"]
                        if pid not in training_ids]

        # Combine
        combined_ids = train_ids_this + ref_ids_this
        n_train = len(train_ids_this)
        n_ref = len(ref_ids_this)
        n_total = len(combined_ids)

        if n_total == 0:
            continue

        print(f"\n  {inc_type}: {n_train:,} training + {n_ref:,} reference "
              f"= {n_total:,} total", flush=True)

        # Build combined vector matrix
        combined_vectors = np.zeros((n_total, N_KMERS), dtype=np.float32)
        for i, pid in enumerate(combined_ids):
            if pid in train_id_to_vec:
                combined_vectors[i] = train_id_to_vec[pid].astype(np.float32)
            elif pid in ref_id_to_idx:
                combined_vectors[i] = ref_vectors[ref_id_to_idx[pid]]

        # Cluster (member keys = plasmid IDs, used for deterministic cluster numbering;
        # inc_type anchors IDs to the frozen training-only run when available)
        assignments = cluster_inc_group(combined_vectors, combined_ids, max_group_size, inc_type=inc_type)

        # Report cluster counts
        for bname in PLIN_THRESHOLDS:
            n_clusters = len(set(assignments[bname]))
            print(f"    Bin {bname}: {n_clusters:>6,} clusters", flush=True)

        # Build pLIN codes and results
        bin_labels = list(PLIN_THRESHOLDS.keys())
        for i, pid in enumerate(combined_ids):
            code_parts = [str(assignments[b][i]) for b in bin_labels]
            plin_code = ".".join(code_parts)

            is_training = pid in training_ids
            source = "training" if is_training else "reference"

            # Get metadata
            if is_training:
                rec = train_id_to_rec[pid]
                length = rec["length"]
                confidence = 1.0
                is_novel = False
                is_multiple = False
                inc_secondary = ""
            else:
                row = ref_rows_this[ref_rows_this["plasmid_id"] == pid].iloc[0]
                length = int(row["length_bp"])
                confidence = float(row["inc_confidence"])
                is_novel = False
                is_multiple = bool(row.get("is_multiple", False))
                inc_secondary = str(row.get("inc_secondary", ""))

            all_results.append({
                "plasmid_id": pid,
                "inc_type": inc_type,
                "resolved_inc_type": resolve_inc_type(inc_type),
                "inc_confidence": confidence,
                "length_bp": length,
                "pLIN": plin_code,
                "bin_A": int(assignments["A"][i]),
                "bin_B": int(assignments["B"][i]),
                "bin_C": int(assignments["C"][i]),
                "bin_D": int(assignments["D"][i]),
                "bin_E": int(assignments["E"][i]),
                "bin_F": int(assignments["F"][i]),
                "source": source,
                "is_multiple": is_multiple,
                "inc_secondary": inc_secondary,
                "is_novel": is_novel,
            })

    elapsed = time.time() - t0

    # Add Unknown/Novel sequences (no pLIN code)
    novel_ref = ref_classifications[ref_classifications["is_novel"]]
    for _, row in novel_ref.iterrows():
        pid = row["plasmid_id"]
        if pid in training_ids:
            continue  # Training sequences always have known Inc type
        all_results.append({
            "plasmid_id": pid,
            "inc_type": "Unknown",
            "resolved_inc_type": "Unknown",
            "inc_confidence": float(row["inc_confidence"]),
            "length_bp": int(row["length_bp"]),
            "pLIN": "NA",
            "bin_A": "NA", "bin_B": "NA", "bin_C": "NA",
            "bin_D": "NA", "bin_E": "NA", "bin_F": "NA",
            "source": "reference",
            "is_multiple": bool(row.get("is_multiple", False)),
            "inc_secondary": "",
            "is_novel": True,
        })

    print(f"\n  Clustering complete in {elapsed/60:.1f} min", flush=True)
    print(f"  Total results: {len(all_results):,}", flush=True)

    return pd.DataFrame(all_results)


# ── Phase 4: Output ──────────────────────────────────────────────────────────

def training_fingerprint():
    """Fingerprint of the training FASTA set, comparable with assign_pLIN.py's."""
    fasta_paths = sorted(glob.glob(os.path.join(TRAINING_DIR, "*", "fastas", "*.fasta")))
    hasher = hashlib.sha256()
    for p in fasta_paths:
        hasher.update(os.path.relpath(p, TRAINING_DIR).encode())
        hasher.update(str(os.path.getsize(p)).encode())
    return hasher.hexdigest()[:12]


def save_results(df):
    """Save final pLIN assignments and print summary."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"\n  Output saved: {OUTPUT_FILE}", flush=True)
    print(f"  Rows: {len(df):,}", flush=True)

    stamp = run_version_stamp(training_fingerprint=training_fingerprint())
    stamp_path = OUTPUT_FILE + ".run_info.json"
    with open(stamp_path, "w") as fh:
        json.dump(stamp, fh, indent=2)
    print(f"  Run info saved: {stamp_path}", flush=True)
    print(f"  Run: {stamp['run_timestamp_utc']}  commit={stamp['git_commit']}  "
          f"training_fingerprint={stamp['training_input_fingerprint']}", flush=True)

    # Summary
    print("\n" + "=" * 70, flush=True)
    print("SUMMARY", flush=True)
    print("=" * 70, flush=True)

    classified = df[df["pLIN"] != "NA"]
    unknown = df[df["pLIN"] == "NA"]
    training = df[df["source"] == "training"]
    reference = df[df["source"] == "reference"]

    print(f"  Total sequences:      {len(df):,}", flush=True)
    print(f"    Training:           {len(training):,}", flush=True)
    print(f"    Reference:          {len(reference):,}", flush=True)
    print(f"  Classified (pLIN):    {len(classified):,} "
          f"({len(classified)/len(df)*100:.1f}%)", flush=True)
    print(f"  Unknown/Novel:        {len(unknown):,} "
          f"({len(unknown)/len(df)*100:.1f}%)", flush=True)

    if len(classified) > 0:
        print(f"  Unique pLIN codes:    {classified['pLIN'].nunique():,}", flush=True)
        print(f"  Inc groups:           {classified['inc_type'].nunique()}", flush=True)

    print(f"\n  Inc type distribution:", flush=True)
    for inc, count in df["inc_type"].value_counts().items():
        n_plin = 0
        if inc != "Unknown":
            n_plin = classified[classified["inc_type"] == inc]["pLIN"].nunique()
        print(f"    {inc:<12s} {count:>6,} sequences, "
              f"{n_plin:>6,} unique pLIN codes", flush=True)

    # Cluster stats for classified sequences
    if len(classified) > 0:
        print(f"\n  Cluster counts by level:", flush=True)
        for b in ["bin_A", "bin_B", "bin_C", "bin_D", "bin_E", "bin_F"]:
            n_clusters = classified[b].nunique()
            print(f"    {b}: {n_clusters:,} clusters", flush=True)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Assign pLIN codes to reference sequences")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from checkpoints if available")
    parser.add_argument("--max-group-size", type=int, default=25000,
                        help="Maximum group size before chunked clustering (default: 25000)")
    args = parser.parse_args()

    print("=" * 70, flush=True)
    print("pLIN Reference Assignment Pipeline", flush=True)
    print(f"  Classifier: {CLASSIFIER_PATH}", flush=True)
    print(f"  Reference:  {SEQUENCES_FASTA if os.path.exists(SEQUENCES_FASTA) else REFERENCE_DIR}", flush=True)
    print(f"  Output:     {OUTPUT_FILE}", flush=True)
    print(f"  Resume:     {args.resume}", flush=True)
    print("=" * 70, flush=True)

    # ── Phase 1: Load sequences & compute 4-mer vectors ──
    print("\n" + "=" * 70, flush=True)
    print("PHASE 1: Load Sequences & Compute 4-mer Vectors", flush=True)
    print("=" * 70, flush=True)

    print("\n[1a] Loading reference sequences ...", flush=True)
    ref_records = load_reference_sequences()

    print("\n[1b] Loading training sequences ...", flush=True)
    training_records = load_training_sequences()

    print("\n[1c] Computing reference 4-mer vectors ...", flush=True)
    ref_vectors, ref_ids, ref_lengths = compute_vectors(
        ref_records, VECTORS_CHECKPOINT, resume=args.resume)

    print("\n[1d] Computing training 4-mer vectors ...", flush=True)
    training_vectors = np.zeros((len(training_records), N_KMERS), dtype=np.float32)
    for i, rec in enumerate(training_records):
        training_vectors[i] = kmer_vector(rec["sequence"]).astype(np.float32)
    print(f"  {len(training_records):,} training vectors computed", flush=True)

    # Free sequence strings to save memory
    for rec in ref_records:
        del rec["sequence"]
    for rec in training_records:
        del rec["sequence"]

    # ── Phase 2: Classify Inc types ──
    print("\n" + "=" * 70, flush=True)
    print("PHASE 2: Classify Inc Types (KNN)", flush=True)
    print("=" * 70, flush=True)

    ref_classifications = classify_inc_types(
        ref_vectors, ref_ids, ref_lengths, INC_CHECKPOINT, resume=args.resume)

    # ── Phase 3: Per-Inc-group clustering ──
    results_df = phase3_clustering(
        ref_classifications, ref_vectors, ref_ids,
        training_records, training_vectors,
        max_group_size=args.max_group_size)

    # ── Phase 4: Output ──
    print("\n" + "=" * 70, flush=True)
    print("PHASE 4: Save Results", flush=True)
    print("=" * 70, flush=True)

    save_results(results_df)

    print("\n" + "=" * 70, flush=True)
    print("Pipeline complete!", flush=True)
    print("=" * 70, flush=True)


if __name__ == "__main__":
    main()
