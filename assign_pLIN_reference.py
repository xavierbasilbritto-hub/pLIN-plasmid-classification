#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Assign pLIN codes to ALL reference plasmid sequences.

Pipeline:
  Phase 1: Compute 4-mer frequency vectors for all reference sequences
  Phase 2: Classify Inc types using the pre-trained KNN classifier
            (reported alongside the code; the code itself does not use it)
  Phase 3: Extend the frozen training founder tree (plin_founder.py) with
            every reference plasmid, in accession order → pLIN codes
  Phase 4: Output combined table

Training codes are copied unchanged from output/pLIN_assignments.tsv and the
founder tree built with them (data/plin_founder_tree_training.npz) is only
ever extended, so no existing code changes when the database grows. Run
assign_pLIN.py first. The extended tree is saved as
data/plin_founder_tree_reference.npz and is what the app's query mode uses.

Usage:
  python assign_pLIN_reference.py [--resume]
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
from sklearn.neighbors import KNeighborsClassifier
from Bio import SeqIO

from plin_founder import PLIN_THRESHOLDS, FounderTree, kmer_vector

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
TRAINING_TREE_FILE = os.path.join(BASE_DIR, "data", "plin_founder_tree_training.npz")
REFERENCE_TREE_FILE = os.path.join(BASE_DIR, "data", "plin_founder_tree_reference.npz")

# Database version file: distinct from PLIN_APP_VERSION (plin_app.py). The
# app version tracks code releases; this tracks database CONTENT, which is
# regenerated on its own schedule (adding sequences, reclassifying, etc.)
# independently of app releases. Lives at the repo root, not under output/,
# so it survives being copied alongside the TSV+classifier into a standalone
# distributable database bundle (see build_database_release.py).
DATABASE_VERSION_PATH = os.path.join(BASE_DIR, "DATABASE_VERSION.json")

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
    """Reproducibility stamp: see assign_pLIN.py::run_version_stamp for rationale."""
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

N_KMERS = 256


# ── Phase 1: Load sequences & compute 4-mer vectors ──────────────────────────

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


# ── Phase 3: Founder assignment ──────────────────────────────────────────────

TRAINING_ASSIGNMENTS_PATH = os.path.join(OUTPUT_DIR, "pLIN_assignments.tsv")


def phase3_assign_codes(ref_classifications, ref_vectors, ref_ids, training_records):
    """Extend the frozen training founder tree with every reference plasmid."""
    print("\n" + "=" * 70, flush=True)
    print("PHASE 3: Founder Assignment (extending the frozen training tree)", flush=True)
    print("=" * 70, flush=True)

    if not (os.path.exists(TRAINING_TREE_FILE) and os.path.exists(TRAINING_ASSIGNMENTS_PATH)):
        sys.exit(f"  Missing {TRAINING_TREE_FILE} or {TRAINING_ASSIGNMENTS_PATH}: run assign_pLIN.py first.")
    tree = FounderTree.from_npz(TRAINING_TREE_FILE)
    training = pd.read_csv(TRAINING_ASSIGNMENTS_PATH, sep="\t")
    bin_cols = [f"bin_{b}" for b in PLIN_THRESHOLDS]

    all_results = []
    train_len = {r["plasmid_id"]: r["length"] for r in training_records}
    for _, row in training.iterrows():
        all_results.append({
            "plasmid_id": row["plasmid_id"], "inc_type": row["inc_type"],
            "resolved_inc_type": resolve_inc_type(row["inc_type"]), "inc_confidence": 1.0,
            "length_bp": int(train_len.get(row["plasmid_id"], row["length_bp"])),
            "pLIN": row["pLIN"], **{c: int(row[c]) for c in bin_cols},
            "source": "training", "is_multiple": False, "inc_secondary": "", "is_novel": False,
        })

    training_ids = set(training["plasmid_id"])
    meta = ref_classifications.set_index("plasmid_id")
    order = sorted(i for i, pid in enumerate(ref_ids) if pid not in training_ids)
    order.sort(key=lambda i: ref_ids[i])
    n_skipped = len(ref_ids) - len(order)
    if n_skipped:
        print(f"  {n_skipped:,} reference sequences are training plasmids "
              f"(training code kept)", flush=True)

    t0 = time.time()
    for n_done, i in enumerate(order, start=1):
        pid = ref_ids[i]
        code, _, _ = tree.assign(ref_vectors[i], key=pid)
        row = meta.loc[pid]
        novel = bool(row["is_novel"])
        inc = "Unknown" if novel else row["inc_type"]
        all_results.append({
            "plasmid_id": pid, "inc_type": inc, "resolved_inc_type": resolve_inc_type(inc),
            "inc_confidence": float(row["inc_confidence"]), "length_bp": int(row["length_bp"]),
            "pLIN": ".".join(map(str, code)), **dict(zip(bin_cols, map(int, code))),
            "source": "reference", "is_multiple": bool(row.get("is_multiple", False)),
            "inc_secondary": "" if novel else str(row.get("inc_secondary", "")),
            "is_novel": novel,
        })
        if n_done % 10000 == 0:
            print(f"    {n_done:,}/{len(order):,} assigned ({time.time() - t0:.0f}s)", flush=True)

    tree.to_npz(REFERENCE_TREE_FILE)
    print(f"  Extended founder tree saved: {REFERENCE_TREE_FILE}", flush=True)
    print(f"  Founder assignment complete in {(time.time() - t0)/60:.1f} min", flush=True)
    return pd.DataFrame(all_results)


def training_fingerprint():
    """Fingerprint of the training FASTA set, comparable with assign_pLIN.py's."""
    fasta_paths = sorted(glob.glob(os.path.join(TRAINING_DIR, "*", "fastas", "*.fasta")))
    hasher = hashlib.sha256()
    for p in fasta_paths:
        hasher.update(os.path.relpath(p, TRAINING_DIR).encode())
        hasher.update(str(os.path.getsize(p)).encode())
    return hasher.hexdigest()[:12]


def write_database_version(df, stamp):
    """Write DATABASE_VERSION.json: a human-facing version identity for the
    reference database, separate from PLIN_APP_VERSION. The version string
    is date-stamped (db-YYYY.MM.DD) rather than incrementing semver, since
    what matters to a user is WHEN this snapshot was built and against what
    input, not a sequential release number: two database builds on the
    same day from different inputs are already disambiguated by
    content_sha256 and training_input_fingerprint below.
    """
    build_date = stamp["run_timestamp_utc"][:10]  # YYYY-MM-DD
    version_string = f"db-{build_date.replace('-', '.')}"

    classified = df[~df["is_novel"].astype(bool)]
    hasher = hashlib.sha256()
    with open(OUTPUT_FILE, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            hasher.update(chunk)

    version_info = {
        "database_version": version_string,
        "build_timestamp_utc": stamp["run_timestamp_utc"],
        "git_commit": stamp["git_commit"],
        "training_input_fingerprint": stamp["training_input_fingerprint"],
        "total_plasmids": int(len(df)),
        "classified_plasmids": int(len(classified)),
        "unique_plin_codes": int(df["pLIN"].nunique()),
        "code_scheme": "founder (plin_founder.py); codes never renumbered across releases",
        "inc_rep_groups": int(df.loc[df["inc_type"] != "Unknown", "inc_type"].nunique()),
        "training_plasmids": int((df["source"] == "training").sum()),
        "reference_plasmids": int((df["source"] == "reference").sum()),
        "content_sha256": hasher.hexdigest(),
        "reference_assignments_file": os.path.basename(OUTPUT_FILE),
        "classifier_file": os.path.basename(CLASSIFIER_PATH),
    }
    with open(DATABASE_VERSION_PATH, "w") as fh:
        json.dump(version_info, fh, indent=2)
    print(f"  Database version saved: {DATABASE_VERSION_PATH}", flush=True)
    print(f"  Database version: {version_string}  "
          f"({version_info['total_plasmids']:,} plasmids, "
          f"{version_info['unique_plin_codes']:,} unique pLIN codes)", flush=True)
    return version_info


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

    write_database_version(df, stamp)

    # Summary
    print("\n" + "=" * 70, flush=True)
    print("SUMMARY", flush=True)
    print("=" * 70, flush=True)

    classified = df[~df["is_novel"].astype(bool)]
    unknown = df[df["is_novel"].astype(bool)]
    training = df[df["source"] == "training"]
    reference = df[df["source"] == "reference"]

    print(f"  Total sequences:      {len(df):,}", flush=True)
    print(f"    Training:           {len(training):,}", flush=True)
    print(f"    Reference:          {len(reference):,}", flush=True)
    print(f"  Unique pLIN codes (all plasmids): {df['pLIN'].nunique():,}", flush=True)
    print(f"  Inc/Rep classified:   {len(classified):,} "
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

    # ── Phase 3: Founder assignment ──
    results_df = phase3_assign_codes(
        ref_classifications, ref_vectors, ref_ids, training_records)

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
