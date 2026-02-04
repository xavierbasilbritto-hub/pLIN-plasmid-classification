#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
pLIN Assignment Script
Assigns plasmid Life Identification Numbers to all sequences in the training folders.
Uses tetranucleotide (4-mer) composition-based cosine distance + single-linkage clustering.
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
TRAINING_DIR = os.path.join(BASE_DIR, "plasmid_sequences_for_training")

INC_TYPES = {
    "IncFII": os.path.join(TRAINING_DIR, "IncFII", "fastas"),
    "IncN":   os.path.join(TRAINING_DIR, "IncN",   "fastas"),
    "IncX1":  os.path.join(TRAINING_DIR, "IncX1",  "fastas"),
}

# pLIN hierarchical thresholds (cosine distance on 4-mer frequencies)
# Mapped from ANI-based thresholds to composition-distance equivalents
PLIN_THRESHOLDS = {
    "A": 0.150,   # ~85% ANI — Broad plasmid family
    "B": 0.100,   # ~90% ANI — Subfamily
    "C": 0.050,   # ~95% ANI — Cluster
    "D": 0.020,   # ~98% ANI — Subcluster
    "E": 0.010,   # ~99% ANI — Clone
    "F": 0.001,   # ~99.9% ANI — Strain / Outbreak
}

OUTPUT_FILE = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")


# ── Step 1: Load sequences ────────────────────────────────────────────────────

def load_all_sequences():
    """Load all plasmid sequences from all Inc type folders."""
    records = []
    for inc_type, fasta_dir in INC_TYPES.items():
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
    print(f"  Total: {len(records)} plasmid sequences\n")
    return records


# ── Step 2: Compute 4-mer composition vectors ─────────────────────────────────

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
        if (idx + 1) % 1000 == 0:
            print(f"    {idx+1}/{len(records)} done")

    print(f"  Vectors shape: {vectors.shape}\n")
    return vectors


# ── Step 3: Pairwise distance + single-linkage clustering ─────────────────────

def assign_plin_codes(records, vectors):
    """
    Cluster plasmids at each pLIN threshold using single-linkage
    and build hierarchical pLIN codes.
    """
    print("  Computing pairwise cosine distances ...")
    dist_condensed = pdist(vectors, metric="cosine")
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

    return plin_codes, cluster_assignments


# ── Step 4: Build results table ───────────────────────────────────────────────

def build_results(records, plin_codes, cluster_assignments):
    """Assemble final results DataFrame."""
    rows = []
    bin_labels = list(PLIN_THRESHOLDS.keys())

    for i, rec in enumerate(records):
        row = {
            "plasmid_id": rec["plasmid_id"],
            "inc_type": rec["inc_type"],
            "length_bp": rec["length"],
            "pLIN": plin_codes[i],
        }
        for b in bin_labels:
            row[f"bin_{b}"] = cluster_assignments[b][i]
        rows.append(row)

    df = pd.DataFrame(rows)
    return df


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("pLIN Assignment — Plasmid Life Identification Numbers")
    print("=" * 70)

    print("\n[1/4] Loading sequences ...")
    records = load_all_sequences()

    print("[2/4] Computing composition vectors ...")
    vectors = compute_kmer_vectors(records, k=4)

    print("[3/4] Assigning pLIN codes ...")
    plin_codes, cluster_assignments = assign_plin_codes(records, vectors)

    print(f"\n[4/4] Building results table ...")
    df = build_results(records, plin_codes, cluster_assignments)

    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    df.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"\n  Results saved to: {OUTPUT_FILE}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Total plasmids:  {len(df)}")
    print(f"  Inc types:       {df['inc_type'].value_counts().to_dict()}")
    print(f"  Unique pLIN codes: {df['pLIN'].nunique()}")
    print()

    # Show per-Inc-type pLIN distribution
    for inc in sorted(df["inc_type"].unique()):
        sub = df[df["inc_type"] == inc]
        print(f"  {inc}: {len(sub)} plasmids → {sub['pLIN'].nunique()} unique pLIN codes")

    print()

    # Show sample assignments
    print("Sample pLIN assignments (first 20):")
    print("-" * 70)
    print(df[["plasmid_id", "inc_type", "length_bp", "pLIN"]].head(20).to_string(index=False))
    print()

    # Show most common pLIN codes
    print("Top 20 most common pLIN codes:")
    print("-" * 70)
    top = df["pLIN"].value_counts().head(20)
    for code, count in top.items():
        inc_dist = df[df["pLIN"] == code]["inc_type"].value_counts().to_dict()
        print(f"  {code:<30s}  n={count:>4d}  {inc_dist}")

    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()
