#!/usr/bin/env python3
"""
pLIN Test Script
Assigns pLIN codes to test plasmid sequences in test_plasmids/ folder.
Uses the same tetranucleotide composition + single-linkage clustering as assign_pLIN.py.
"""

import os
import glob
import numpy as np
import pandas as pd
from itertools import product as iter_product
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
from Bio import SeqIO

# ── Configuration ──────────────────────────────────────────────────────────────

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.join(BASE_DIR, "test_plasmids")

# Auto-detect Inc-group subfolders that contain FASTA files
TEST_INC_TYPES = {}
if os.path.isdir(TEST_DIR):
    for entry in sorted(os.listdir(TEST_DIR)):
        subdir = os.path.join(TEST_DIR, entry)
        if os.path.isdir(subdir):
            fastas = glob.glob(os.path.join(subdir, "*.fasta"))
            if fastas:
                TEST_INC_TYPES[entry] = subdir

# pLIN hierarchical thresholds (same as assign_pLIN.py)
PLIN_THRESHOLDS = {
    "A": 0.150,   # ~85% ANI — Broad plasmid family
    "B": 0.100,   # ~90% ANI — Subfamily
    "C": 0.050,   # ~95% ANI — Cluster
    "D": 0.020,   # ~98% ANI — Subcluster
    "E": 0.010,   # ~99% ANI — Clone
    "F": 0.001,   # ~99.9% ANI — Strain / Outbreak
}

OUTPUT_DIR = os.path.join(BASE_DIR, "output", "test")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "test_pLIN_assignments.tsv")


# ── Step 1: Load test sequences ──────────────────────────────────────────────

def load_test_sequences():
    """Load all plasmid sequences from test_plasmids/ subfolders."""
    records = []
    for inc_type, fasta_dir in TEST_INC_TYPES.items():
        fasta_files = sorted(glob.glob(os.path.join(fasta_dir, "*.fasta")))
        count = 0
        for fasta_file in fasta_files:
            for rec in SeqIO.parse(fasta_file, "fasta"):
                records.append({
                    "plasmid_id": rec.id,
                    "inc_type": inc_type,
                    "sequence": str(rec.seq),
                    "length": len(rec.seq),
                    "source_file": os.path.basename(fasta_file),
                })
                count += 1
        print(f"  Loaded {count:>5} sequences from {inc_type}")
    print(f"  Total: {len(records)} test plasmid sequences\n")
    return records


# ── Step 2: Compute 4-mer composition vectors ────────────────────────────────

def compute_kmer_vectors(records, k=4):
    """Compute normalised tetranucleotide frequency vectors."""
    bases = "ACGT"
    all_kmers = ["".join(p) for p in iter_product(bases, repeat=k)]

    print(f"  Computing {k}-mer frequency vectors for {len(records)} plasmids ...")
    vectors = np.zeros((len(records), len(all_kmers)), dtype=np.float64)

    for idx, rec in enumerate(records):
        seq = rec["sequence"].upper()
        total = max(len(seq) - k + 1, 1)
        for ki, kmer in enumerate(all_kmers):
            vectors[idx, ki] = seq.count(kmer) / total
        if (idx + 1) % 100 == 0:
            print(f"    {idx+1}/{len(records)} done")

    print(f"  Vectors shape: {vectors.shape}\n")
    return vectors


# ── Step 3: Pairwise distance + single-linkage clustering ────────────────────

def assign_plin_codes(records, vectors):
    """Cluster test plasmids at each pLIN threshold using single-linkage."""
    print("  Computing pairwise cosine distances ...")
    dist_condensed = pdist(vectors, metric="cosine")
    dist_matrix = squareform(dist_condensed)
    print(f"  {len(dist_condensed)} pairwise distances computed\n")

    Z = linkage(dist_condensed, method="single")

    bin_labels = list(PLIN_THRESHOLDS.keys())
    bin_thresholds = list(PLIN_THRESHOLDS.values())

    cluster_assignments = {}
    for bname, thresh in zip(bin_labels, bin_thresholds):
        clusters = fcluster(Z, t=thresh, criterion="distance")
        cluster_assignments[bname] = clusters
        n_clusters = len(set(clusters))
        print(f"  Bin {bname} (d ≤ {thresh:.3f}): {n_clusters:>5} clusters")

    # Build pLIN codes
    n = len(records)
    plin_codes = []
    for i in range(n):
        code_parts = [str(cluster_assignments[b][i]) for b in bin_labels]
        plin_codes.append(".".join(code_parts))

    return plin_codes, cluster_assignments, dist_matrix


# ── Step 4: Build results table ──────────────────────────────────────────────

def build_results(records, plin_codes, cluster_assignments):
    """Assemble final results DataFrame."""
    rows = []
    bin_labels = list(PLIN_THRESHOLDS.keys())

    for i, rec in enumerate(records):
        row = {
            "plasmid_id": rec["plasmid_id"],
            "inc_type": rec["inc_type"],
            "source_file": rec["source_file"],
            "length_bp": rec["length"],
            "pLIN": plin_codes[i],
        }
        for b in bin_labels:
            row[f"bin_{b}"] = cluster_assignments[b][i]
        rows.append(row)

    df = pd.DataFrame(rows)
    return df


# ── Step 5: Distance summary ────────────────────────────────────────────────

def print_distance_summary(records, dist_matrix):
    """Print pairwise distance statistics."""
    n = len(records)
    upper_tri = dist_matrix[np.triu_indices(n, k=1)]

    print(f"\n  Pairwise cosine distance statistics:")
    print(f"    Min:    {upper_tri.min():.6f}")
    print(f"    Max:    {upper_tri.max():.6f}")
    print(f"    Mean:   {upper_tri.mean():.6f}")
    print(f"    Median: {np.median(upper_tri):.6f}")
    print(f"    Std:    {upper_tri.std():.6f}")

    # How many pairs fall within each threshold
    print(f"\n  Pairs within each pLIN threshold:")
    total_pairs = len(upper_tri)
    for bname, thresh in PLIN_THRESHOLDS.items():
        count = np.sum(upper_tri <= thresh)
        pct = 100.0 * count / total_pairs
        print(f"    Bin {bname} (d ≤ {thresh:.3f}): {count:>6} / {total_pairs} pairs ({pct:.1f}%)")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("pLIN TEST — Assigning pLIN codes to test plasmids")
    print("=" * 70)

    if not TEST_INC_TYPES:
        print("\n  ERROR: No FASTA files found in test_plasmids/ subfolders.")
        print("  Place .fasta files in test_plasmids/<IncGroup>/ directories.")
        return

    print(f"\n  Test directories with FASTA files:")
    for inc, path in TEST_INC_TYPES.items():
        n = len(glob.glob(os.path.join(path, "*.fasta")))
        print(f"    {inc}: {n} FASTA files")

    print(f"\n[1/4] Loading test sequences ...")
    records = load_test_sequences()

    if len(records) < 2:
        print("  ERROR: Need at least 2 sequences for clustering.")
        return

    print("[2/4] Computing composition vectors ...")
    vectors = compute_kmer_vectors(records, k=4)

    print("[3/4] Assigning pLIN codes ...")
    plin_codes, cluster_assignments, dist_matrix = assign_plin_codes(records, vectors)

    print(f"\n[4/4] Building results table ...")
    df = build_results(records, plin_codes, cluster_assignments)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"\n  Results saved to: {OUTPUT_FILE}")

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Total test plasmids:   {len(df)}")
    print(f"  Inc types:             {df['inc_type'].value_counts().to_dict()}")
    print(f"  Unique pLIN codes:     {df['pLIN'].nunique()}")

    for inc in sorted(df["inc_type"].unique()):
        sub = df[df["inc_type"] == inc]
        print(f"\n  {inc}: {len(sub)} plasmids → {sub['pLIN'].nunique()} unique pLIN codes")

    # Distance statistics
    print_distance_summary(records, dist_matrix)

    # Full results table
    print("\n" + "=" * 70)
    print("ALL TEST pLIN ASSIGNMENTS")
    print("=" * 70)
    print(df[["plasmid_id", "inc_type", "source_file", "length_bp", "pLIN"]].to_string(index=False))

    # pLIN code distribution
    print("\n" + "-" * 70)
    print("pLIN code distribution:")
    print("-" * 70)
    for code, count in df["pLIN"].value_counts().items():
        members = df[df["pLIN"] == code]["source_file"].tolist()
        print(f"  {code:<30s}  n={count:>2d}  {members}")

    # Nearest-neighbour pairs
    print("\n" + "-" * 70)
    print("Closest plasmid pairs (top 10):")
    print("-" * 70)
    n = len(records)
    pairs = []
    for i in range(n):
        for j in range(i+1, n):
            pairs.append((dist_matrix[i, j], records[i]["source_file"], records[j]["source_file"]))
    pairs.sort()
    for dist, a, b in pairs[:10]:
        print(f"  d={dist:.6f}  {a} ↔ {b}")

    print("\n" + "=" * 70)
    print("Test complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
