#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Assign pLIN codes to ALL reference plasmid sequences (~72,556).

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
import time
import argparse
import numpy as np
import pandas as pd
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
MULTI_INC_THRESHOLD = 0.25

# 4-mer setup
K = 4
BASES = "ACGT"
KMERS = ["".join(p) for p in iter_product(BASES, repeat=K)]
KMER_INDEX = {km: i for i, km in enumerate(KMERS)}
N_KMERS = len(KMERS)  # 256


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
                    "inc_secondary": "",
                })
    print(f"  Training sequences: {len(records):,} across "
          f"{len(set(r['inc_type'] for r in records))} Inc groups", flush=True)
    return records


def compute_vectors(records, checkpoint_path, resume=False):
    """Compute 4-mer vectors for all records, with checkpoint support."""
    if resume and os.path.exists(checkpoint_path):
        print(f"  Resuming: loading vectors from {checkpoint_path} ...", flush=True)
        data = np.load(checkpoint_path)
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

def cluster_inc_group(vectors, max_group_size=25000):
    """Run single-linkage clustering on vectors within one Inc group."""
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

    # Cut at each threshold
    assignments = {}
    for bname, thresh in PLIN_THRESHOLDS.items():
        clusters = fcluster(Z, t=thresh, criterion="distance")
        assignments[bname] = clusters

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

        # Cluster
        assignments = cluster_inc_group(combined_vectors, max_group_size)

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
                inc_secondary = ""
            else:
                row = ref_rows_this[ref_rows_this["plasmid_id"] == pid].iloc[0]
                length = int(row["length_bp"])
                confidence = float(row["inc_confidence"])
                is_novel = False
                inc_secondary = str(row.get("inc_secondary", ""))

            all_results.append({
                "plasmid_id": pid,
                "inc_type": inc_type,
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
            "inc_confidence": float(row["inc_confidence"]),
            "length_bp": int(row["length_bp"]),
            "pLIN": "NA",
            "bin_A": "NA", "bin_B": "NA", "bin_C": "NA",
            "bin_D": "NA", "bin_E": "NA", "bin_F": "NA",
            "source": "reference",
            "inc_secondary": "",
            "is_novel": True,
        })

    print(f"\n  Clustering complete in {elapsed/60:.1f} min", flush=True)
    print(f"  Total results: {len(all_results):,}", flush=True)

    return pd.DataFrame(all_results)


# ── Phase 4: Output ──────────────────────────────────────────────────────────

def save_results(df):
    """Save final pLIN assignments and print summary."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"\n  Output saved: {OUTPUT_FILE}", flush=True)
    print(f"  Rows: {len(df):,}", flush=True)

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
