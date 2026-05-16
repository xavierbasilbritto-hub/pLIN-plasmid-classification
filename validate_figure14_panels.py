#!/usr/bin/env python3
"""Compute real data for Figure 14 Panels C and D.

Panel C: Assembly completeness scores for training plasmids
Panel D: Leave-one-out nearest-neighbour distances (real NN distance distribution)
"""

import os
import sys
import json
import glob
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "output", "figure14_validation")
os.makedirs(OUT_DIR, exist_ok=True)

CLASSIFIER_PATH = os.path.join(BASE_DIR, "data", "inc_classifier.npz")
TRAINING_DIR = os.path.join(BASE_DIR, "plasmid_sequences_for_training")


def compute_assembly_completeness():
    """Run assess_assembly_completeness on all training FASTAs."""
    import re
    from Bio import SeqIO

    print("=" * 60)
    print("Panel C: Computing assembly completeness for training plasmids")
    print("=" * 60)

    results = []
    fasta_dirs = sorted(glob.glob(os.path.join(TRAINING_DIR, "*/fastas")))
    total = 0

    for fasta_dir in fasta_dirs:
        inc_group = os.path.basename(os.path.dirname(fasta_dir))
        fastas = sorted(glob.glob(os.path.join(fasta_dir, "*.fasta")) +
                       glob.glob(os.path.join(fasta_dir, "*.fa")) +
                       glob.glob(os.path.join(fasta_dir, "*.fna")))

        for fpath in fastas:
            try:
                rec = next(SeqIO.parse(fpath, "fasta"))
                seq = str(rec.seq).upper()
                total_len = len(seq)
                pid = rec.id

                if total_len == 0:
                    continue

                # 1. Contig count
                contigs = re.split(r'N{10,}', seq)
                contigs = [c for c in contigs if len(c) > 0]
                n_contigs = max(len(contigs), 1)
                n_gaps = max(n_contigs - 1, 0)

                # 2. N50 ratio
                contig_lengths = sorted([len(c) for c in contigs], reverse=True)
                cumsum = 0
                n50 = contig_lengths[0]
                half = total_len / 2
                for cl in contig_lengths:
                    cumsum += cl
                    if cumsum >= half:
                        n50 = cl
                        break
                n50_ratio = n50 / total_len if total_len > 0 else 0

                # 3. Circular topology signal
                circular_signal = False
                check_len = min(500, total_len // 4)
                if check_len >= 50:
                    head = seq[:check_len]
                    tail = seq[-check_len:]
                    matches = sum(1 for a, b in zip(head, tail) if a == b)
                    identity = matches / check_len
                    circular_signal = identity >= 0.90

                # 4. No Prodigal data available — neutral score
                coding_density = 0.0

                # 5. Composite score
                score = 0
                score += 40 if n_contigs == 1 else max(0, 40 - (n_contigs - 1) * 10)
                score += 20 if n50_ratio > 0.9 else int(20 * n50_ratio)
                score += 20 if circular_signal else 0
                score += 5  # neutral (no Prodigal)
                score += 10 if n_gaps == 0 else max(0, 10 - n_gaps * 2)
                score = min(100, max(0, score))

                if score >= 80:
                    status = "COMPLETE"
                elif score >= 60:
                    status = "NEAR-COMPLETE"
                elif score >= 40:
                    status = "FRAGMENTED"
                else:
                    status = "POOR"

                results.append({
                    "plasmid_id": pid,
                    "inc_group": inc_group,
                    "total_length": total_len,
                    "n_contigs": n_contigs,
                    "n50_ratio": round(n50_ratio, 3),
                    "circular_signal": circular_signal,
                    "n_gaps": n_gaps,
                    "completeness_score": score,
                    "completeness_status": status,
                })
                total += 1
                if total % 500 == 0:
                    print(f"  Processed {total} plasmids ...")

            except Exception as e:
                pass  # skip problematic files

    df = pd.DataFrame(results)
    out_path = os.path.join(OUT_DIR, "assembly_completeness.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\n  Total: {len(df)} plasmids assessed")
    print(f"  Saved: {out_path}")

    # Summary
    for status in ["COMPLETE", "NEAR-COMPLETE", "FRAGMENTED", "POOR"]:
        n = (df["completeness_status"] == status).sum()
        pct = n / len(df) * 100
        print(f"    {status}: {n} ({pct:.1f}%)")

    print(f"  Score median: {df['completeness_score'].median():.0f}")
    print(f"  Score mean: {df['completeness_score'].mean():.1f}")
    print(f"  Circular signal: {df['circular_signal'].sum()} ({df['circular_signal'].mean()*100:.1f}%)")
    print(f"  Single contig: {(df['n_contigs'] == 1).sum()} ({(df['n_contigs'] == 1).mean()*100:.1f}%)")

    return df


def compute_nn_distances():
    """Compute leave-one-out NN distances for all training samples."""
    print("\n" + "=" * 60)
    print("Panel D: Computing leave-one-out NN distances")
    print("=" * 60)

    data = np.load(CLASSIFIER_PATH, allow_pickle=True)
    X = data["X"].astype(np.float64)
    y = data["y"]
    group_names = list(data["group_names"])

    n_samples = len(X)
    print(f"  Training samples: {n_samples}")
    print(f"  Groups: {len(group_names)}")

    # Compute pairwise distances within each group (leave-one-out NN)
    all_results = []

    for gi, gname in enumerate(group_names):
        mask = y == gi
        X_group = X[mask]
        n_group = len(X_group)

        if n_group < 2:
            continue

        # Compute pairwise cosine distances within group
        dists = cdist(X_group, X_group, metric="cosine")
        # Set diagonal to inf so we don't pick self
        np.fill_diagonal(dists, np.inf)

        # Leave-one-out NN distance for each sample
        nn_dists = dists.min(axis=1)

        for j, d in enumerate(nn_dists):
            # Novelty flag
            if d > 0.050:
                novelty = "Potentially novel lineage"
                indicator = "RED"
            elif d > 0.020:
                novelty = "Divergent"
                indicator = "YELLOW"
            else:
                novelty = "None"
                indicator = "GREEN"

            all_results.append({
                "inc_group": gname,
                "nn_distance": round(float(d), 6),
                "novelty_flag": novelty,
                "coverage_indicator": indicator,
            })

        print(f"  {gname}: n={n_group}, median NN dist={np.median(nn_dists):.4f}, "
              f"max={np.max(nn_dists):.4f}")

    df = pd.DataFrame(all_results)
    out_path = os.path.join(OUT_DIR, "nn_distances.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\n  Total: {len(df)} samples")
    print(f"  Saved: {out_path}")

    # Summary
    for indicator in ["GREEN", "YELLOW", "RED"]:
        n = (df["coverage_indicator"] == indicator).sum()
        pct = n / len(df) * 100
        print(f"    {indicator}: {n} ({pct:.1f}%)")

    print(f"  Overall median NN dist: {df['nn_distance'].median():.4f}")
    print(f"  Overall mean NN dist: {df['nn_distance'].mean():.4f}")
    print(f"  p95 NN dist: {df['nn_distance'].quantile(0.95):.4f}")
    print(f"  p99 NN dist: {df['nn_distance'].quantile(0.99):.4f}")

    # Save summary JSON
    summary = {
        "n_samples": len(df),
        "n_groups": len(group_names),
        "median_nn_dist": round(float(df["nn_distance"].median()), 4),
        "mean_nn_dist": round(float(df["nn_distance"].mean()), 4),
        "p95_nn_dist": round(float(df["nn_distance"].quantile(0.95)), 4),
        "p99_nn_dist": round(float(df["nn_distance"].quantile(0.99)), 4),
        "n_green": int((df["coverage_indicator"] == "GREEN").sum()),
        "n_yellow": int((df["coverage_indicator"] == "YELLOW").sum()),
        "n_red": int((df["coverage_indicator"] == "RED").sum()),
        "pct_green": round(float((df["coverage_indicator"] == "GREEN").mean() * 100), 1),
        "pct_yellow": round(float((df["coverage_indicator"] == "YELLOW").mean() * 100), 1),
        "pct_red": round(float((df["coverage_indicator"] == "RED").mean() * 100), 1),
    }

    summary_path = os.path.join(OUT_DIR, "validation_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  Saved: {summary_path}")

    return df


if __name__ == "__main__":
    comp_df = compute_assembly_completeness()
    nn_df = compute_nn_distances()

    # Save combined summary
    comp_summary = {
        "n_plasmids": len(comp_df),
        "score_median": float(comp_df["completeness_score"].median()),
        "score_mean": round(float(comp_df["completeness_score"].mean()), 1),
        "n_complete": int((comp_df["completeness_status"] == "COMPLETE").sum()),
        "n_near_complete": int((comp_df["completeness_status"] == "NEAR-COMPLETE").sum()),
        "n_fragmented": int((comp_df["completeness_status"] == "FRAGMENTED").sum()),
        "n_poor": int((comp_df["completeness_status"] == "POOR").sum()),
        "pct_complete": round(float((comp_df["completeness_status"] == "COMPLETE").mean() * 100), 1),
        "pct_circular": round(float(comp_df["circular_signal"].mean() * 100), 1),
        "pct_single_contig": round(float((comp_df["n_contigs"] == 1).mean() * 100), 1),
    }

    combined = {
        "completeness": comp_summary,
        "nn_distances": json.load(open(os.path.join(OUT_DIR, "validation_summary.json"))),
    }
    combined_path = os.path.join(OUT_DIR, "figure14_summary.json")
    with open(combined_path, "w") as f:
        json.dump(combined, f, indent=2)
    print(f"\n  Combined summary: {combined_path}")
    print("\nDone!")
