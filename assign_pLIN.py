#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
pLIN Assignment Script
Assigns plasmid Lineage Identification Numbers to all sequences in the training folders.
Uses tetranucleotide (4-mer) composition-based cosine distance + incremental
founder assignment (plin_founder.py): plasmids are processed in accession
order and, at each level, join the nearest cluster founder within that level's
threshold or found a new cluster. Codes never change when plasmids are added,
and clusters cannot chain (see plin_founder.py and benchmark_plin_schemes.py).

The founder tree is saved to data/plin_founder_tree_training.npz; the
reference pipeline (assign_pLIN_reference.py) extends that tree rather than
re-clustering, so training codes are preserved in every release.
"""

import os
import glob
import hashlib
import subprocess
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from Bio import SeqIO

from plin_founder import PLIN_THRESHOLDS, FounderTree, canonical_order, kmer_vector

# ── Configuration ──────────────────────────────────────────────────────────────

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TRAINING_DIR = os.path.join(BASE_DIR, "plasmid_sequences_for_training")

def _discover_inc_types():
    """Dynamically discover all Inc type training folders."""
    inc_types = {}
    for inc_dir in sorted(glob.glob(os.path.join(TRAINING_DIR, "*", "fastas"))):
        inc_name = os.path.basename(os.path.dirname(inc_dir))
        if glob.glob(os.path.join(inc_dir, "*.fasta")):
            inc_types[inc_name] = inc_dir
    return inc_types

INC_TYPES = _discover_inc_types()

OUTPUT_FILE = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")
FOUNDER_TREE_FILE = os.path.join(BASE_DIR, "data", "plin_founder_tree_training.npz")

# Training-folder provenance labels retained as-is in `inc_type` for
# traceability. `resolved_inc_type` maps them onto current PlasmidFinder
# nomenclature, since some legacy labels (e.g. IncAC2) are no longer used
# by PlasmidFinder: IncA and IncC are now typed as separate, non-compatible
# replicon markers (Carattoli, Int. J. Med. Microbiol. 2013; PMID: 29486211),
# so a plasmid trained under the legacy "IncAC2" folder resolves to "IncA/IncC"
# to make clear it predates that split rather than representing a current
# PlasmidFinder marker.
RESOLVED_INC_TYPE = {
    "IncAC2": "IncA/IncC",
}


def resolve_inc_type(inc_type):
    """Map a training-provenance Inc label onto current PlasmidFinder nomenclature."""
    return RESOLVED_INC_TYPE.get(inc_type, inc_type)


def run_version_stamp():
    """Build a reproducibility stamp for this pipeline run.

    Re-running the clustering pipeline (e.g. after adding training data,
    or after any code change) can change pLIN codes even for previously
    assigned plasmids. Stamping every output with a run ID and the exact
    input-file fingerprint makes it possible to tell, at a glance, whether
    two tables/figures were generated from the same run.
    """
    try:
        git_hash = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=BASE_DIR,
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        git_hash = "unknown"

    fasta_paths = sorted(glob.glob(os.path.join(TRAINING_DIR, "*", "fastas", "*.fasta")))
    hasher = hashlib.sha256()
    for p in fasta_paths:
        hasher.update(os.path.relpath(p, TRAINING_DIR).encode())
        hasher.update(str(os.path.getsize(p)).encode())
    input_fingerprint = hasher.hexdigest()[:12]

    return {
        "run_timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "git_commit": git_hash,
        "input_fingerprint": input_fingerprint,
        "n_input_fastas": len(fasta_paths),
    }


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
    """Compute tetranucleotide frequency vectors (plin_founder.kmer_vector)."""
    print(f"  Computing {k}-mer frequency vectors for {len(records)} plasmids ...")
    vectors = np.zeros((len(records), 4 ** k), dtype=np.float64)
    for idx, rec in enumerate(records):
        vectors[idx] = kmer_vector(rec["sequence"], k)
        if (idx + 1) % 1000 == 0:
            print(f"    {idx+1}/{len(records)} done")
    print(f"  Vectors shape: {vectors.shape}\n")
    return vectors


# ── Step 3: Founder assignment ────────────────────────────────────────────────

def assign_plin_codes(records, vectors):
    """Assign founder-based pLIN codes in canonical (accession) order.

    A plasmid ID present in two training folders (21 multi-replicon
    E. faecium plasmids) has an identical vector in both, so its second
    occurrence joins the first one's clusters at distance 0.
    """
    bin_labels = list(PLIN_THRESHOLDS.keys())
    keys = [rec["plasmid_id"] for rec in records]
    tree = FounderTree(PLIN_THRESHOLDS)
    codes = [None] * len(records)
    for n_done, i in enumerate(canonical_order(keys), start=1):
        codes[i], _, _ = tree.assign(vectors[i], key=keys[i])
        if n_done % 2000 == 0:
            print(f"    {n_done}/{len(records)} assigned")

    cluster_assignments = {b: np.array([c[j] for c in codes]) for j, b in enumerate(bin_labels)}
    for b, thresh in PLIN_THRESHOLDS.items():
        print(f"  Bin {b} (d ≤ {thresh:.3f}): {len(set(cluster_assignments[b])):>5} clusters")

    os.makedirs(os.path.dirname(FOUNDER_TREE_FILE), exist_ok=True)
    tree.to_npz(FOUNDER_TREE_FILE)
    print(f"  Founder tree saved to: {FOUNDER_TREE_FILE}")

    plin_codes = [".".join(map(str, c)) for c in codes]
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
            "resolved_inc_type": resolve_inc_type(rec["inc_type"]),
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
    print("pLIN Assignment: Plasmid Lineage Identification Numbers")
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

    stamp = run_version_stamp()
    stamp_path = OUTPUT_FILE + ".run_info.json"
    import json
    with open(stamp_path, "w") as fh:
        json.dump(stamp, fh, indent=2)
    print(f"  Run info saved to: {stamp_path}")
    print(f"  Run: {stamp['run_timestamp_utc']}  "
          f"commit={stamp['git_commit']}  "
          f"input_fingerprint={stamp['input_fingerprint']}")

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
