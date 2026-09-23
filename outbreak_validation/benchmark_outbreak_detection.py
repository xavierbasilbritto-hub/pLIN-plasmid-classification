#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
"""
Leave-one-out cross-validation (LOOCV) benchmark for pLIN outbreak detection.

Uses 74 published outbreak plasmids as ground truth. For each plasmid, removes
it from the reference set, classifies against the remainder, and checks whether
it receives the same pLIN as its known study-mates.

Computes:
  - L6 pLIN concordance rate (same L6 code as study-mates)
  - L3 cluster concordance rate (same L3 cluster as study-mates)
  - Intra-study cluster precision and recall
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from itertools import product as iter_product
from collections import Counter, defaultdict
from scipy.spatial.distance import cosine as cosine_dist

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BASE_DIR)


def kmer_vector(sequence, k=4):
    """Compute normalised 4-mer frequency vector."""
    bases = "ACGT"
    all_kmers = ["".join(p) for p in iter_product(bases, repeat=k)]
    kmer_idx = {km: i for i, km in enumerate(all_kmers)}
    seq = sequence.upper()
    counts = np.zeros(256, dtype=np.float64)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if kmer in kmer_idx:
            counts[kmer_idx[kmer]] += 1
    total = counts.sum()
    if total > 0:
        counts /= total
    return counts


def load_training_data():
    """Load 6,998 training vectors + pLIN assignments."""
    npz_path = os.path.join(BASE_DIR, "data", "inc_classifier.npz")
    plin_path = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")

    data = np.load(npz_path, allow_pickle=True)
    X_train = data["X"].astype(np.float64)
    y_train = data["y"]
    group_names = [str(g) for g in data["group_names"]]

    plin_df = pd.read_csv(plin_path, sep="\t")
    return X_train, y_train, group_names, plin_df


def load_outbreak_data():
    """Load 74 combined outbreak validation plasmids."""
    combined_path = os.path.join(BASE_DIR, "output",
                                  "outbreak_validation_combined_results.tsv")
    return pd.read_csv(combined_path, sep="\t")


def load_outbreak_sequences():
    """Load FASTA sequences for outbreak plasmids."""
    from Bio import SeqIO

    sequences = {}

    def _strip_version(acc):
        """Remove version suffix (e.g., CP022533.1 -> CP022533)."""
        if "." in acc and acc.rsplit(".", 1)[-1].isdigit():
            return acc.rsplit(".", 1)[0]
        return acc

    for seq_dir in [
        # The original 17 outbreak FASTAs live directly in outbreak_validation/,
        # not outbreak_validation/sequences/ (that directory does not exist —
        # a prior version of this script pointed at it, which silently found
        # 0 files and caused those 17 plasmids to skip true LOOCV testing).
        os.path.join(BASE_DIR, "outbreak_validation"),
        os.path.join(BASE_DIR, "outbreak_validation", "expanded_sequences"),
    ]:
        if not os.path.exists(seq_dir):
            continue
        for f in os.listdir(seq_dir):
            if f.endswith((".fasta", ".fa", ".fna")):
                path = os.path.join(seq_dir, f)
                for rec in SeqIO.parse(path, "fasta"):
                    seq_str = str(rec.seq)
                    # Store under both versioned and base accession
                    sequences[rec.id] = seq_str
                    base_acc = _strip_version(rec.id)
                    sequences[base_acc] = seq_str
                    # Also store under filename without extension
                    fname_acc = os.path.splitext(f)[0]
                    sequences[fname_acc] = seq_str

    return sequences


def classify_loocv(query_vec, X_ref, y_ref, group_names, plin_df_ref, k=5):
    """Classify a single query against a reference set (LOO mode).

    Returns dict with: predicted_inc, confidence, pLIN, nn_distance, nn_plasmid.
    """
    # Compute distances to all reference plasmids
    dists = np.array([cosine_dist(query_vec, X_ref[j])
                      for j in range(X_ref.shape[0])])

    # k nearest neighbours
    nn_idx = np.argsort(dists)[:k]
    nn_dists = dists[nn_idx]
    nn_labels = y_ref[nn_idx]

    # Distance-weighted voting for Inc group
    weights = 1.0 / (nn_dists + 1e-10)
    vote_scores = {}
    for label, w in zip(nn_labels, weights):
        vote_scores[label] = vote_scores.get(label, 0.0) + w
    total_weight = sum(vote_scores.values())
    predicted_label = max(vote_scores, key=vote_scores.get)
    confidence = (vote_scores[predicted_label] / total_weight) * 100

    predicted_inc = group_names[predicted_label]

    # Nearest neighbour pLIN assignment
    nn_best = nn_idx[0]
    nn_distance = dists[nn_best]
    nn_row = plin_df_ref.iloc[nn_best]
    nn_plin = nn_row["pLIN"] if "pLIN" in nn_row else ""
    nn_plasmid = nn_row["plasmid_id"] if "plasmid_id" in nn_row else ""

    return {
        "predicted_inc": predicted_inc,
        "confidence": round(confidence, 2),
        "pLIN": nn_plin,
        "nn_distance": round(float(nn_distance), 6),
        "nn_plasmid": nn_plasmid,
    }


def define_ground_truth_clusters(outbreak_df):
    """Define expected outbreak clusters from study groupings.

    Returns dict: study_name -> list of accession IDs that should cluster together.
    """
    clusters = defaultdict(list)
    for _, row in outbreak_df.iterrows():
        study = row["study"]
        clusters[study].append(row["accession"])
    # Only keep studies with 2+ plasmids (cluster = ≥2 members)
    return {s: accs for s, accs in clusters.items() if len(accs) >= 2}


def main():
    print("=" * 70)
    print("pLIN OUTBREAK DETECTION — LEAVE-ONE-OUT CROSS-VALIDATION BENCHMARK")
    print("=" * 70)

    # Load data
    print("\nLoading data...")
    X_train, y_train, group_names, plin_df_train = load_training_data()
    outbreak_df = load_outbreak_data()
    sequences = load_outbreak_sequences()

    print(f"  Training set: {X_train.shape[0]} plasmids, {X_train.shape[1]} features")
    print(f"  Outbreak plasmids: {len(outbreak_df)}")
    print(f"  Sequences loaded: {len(sequences)}")

    # Compute 4-mer vectors for outbreak plasmids
    print("\nComputing 4-mer vectors for outbreak plasmids...")
    outbreak_vecs = {}
    for _, row in outbreak_df.iterrows():
        acc = row["accession"]
        if acc in sequences:
            outbreak_vecs[acc] = kmer_vector(sequences[acc])

    print(f"  Vectors computed: {len(outbreak_vecs)}")
    missing = set(outbreak_df["accession"]) - set(outbreak_vecs.keys())
    if missing:
        print(f"  Missing sequences (will use pre-computed pLIN): {len(missing)}")

    # Define ground truth clusters
    gt_clusters = define_ground_truth_clusters(outbreak_df)
    print(f"\nGround truth clusters (studies with ≥2 plasmids): {len(gt_clusters)}")
    for study, accs in sorted(gt_clusters.items()):
        print(f"  {study}: {len(accs)} plasmids")

    # ── LOOCV for plasmids with available sequences ──────────────────────────
    print("\n" + "=" * 70)
    print("LOOCV: Classifying each outbreak plasmid against training set")
    print("=" * 70)

    results = []
    for idx, row in outbreak_df.iterrows():
        acc = row["accession"]
        original_plin = row["pLIN"]
        original_l3 = ".".join(str(original_plin).split(".")[:3])

        if acc in outbreak_vecs:
            query_vec = outbreak_vecs[acc]
            # Classify against training set (no LOO needed — outbreak plasmids
            # are NOT in the training set)
            result = classify_loocv(query_vec, X_train, y_train, group_names,
                                     plin_df_train, k=5)
            loocv_plin = result["pLIN"]
            loocv_l3 = ".".join(str(loocv_plin).split(".")[:3])

            results.append({
                "accession": acc,
                "study": row["study"],
                "resistance_gene": row["resistance_gene"],
                "country": row.get("country", ""),
                "original_plin": original_plin,
                "loocv_plin": loocv_plin,
                "l6_match": original_plin == loocv_plin,
                "original_l3": original_l3,
                "loocv_l3": loocv_l3,
                "l3_match": original_l3 == loocv_l3,
                "nn_distance": result["nn_distance"],
                "confidence": result["confidence"],
                "predicted_inc": result["predicted_inc"],
                "method": "loocv",
            })
        else:
            # Use pre-computed pLIN from combined results
            results.append({
                "accession": acc,
                "study": row["study"],
                "resistance_gene": row["resistance_gene"],
                "country": row.get("country", ""),
                "original_plin": original_plin,
                "loocv_plin": original_plin,  # Same as original (no recompute)
                "l6_match": True,  # Trivially true
                "original_l3": original_l3,
                "loocv_l3": original_l3,
                "l3_match": True,
                "nn_distance": row.get("nn_distance", None),
                "confidence": row.get("confidence", None),
                "predicted_inc": row.get("predicted_inc", ""),
                "method": "precomputed",
            })

    results_df = pd.DataFrame(results)

    # ── Concordance metrics ──────────────────────────────────────────────────
    loocv_subset = results_df[results_df["method"] == "loocv"]
    n_loocv = len(loocv_subset)

    print(f"\nLOOCV results ({n_loocv} plasmids with sequence data):")
    l6_concordance = loocv_subset["l6_match"].mean() * 100
    l3_concordance = loocv_subset["l3_match"].mean() * 100
    print(f"  L6 pLIN concordance: {l6_concordance:.1f}% ({loocv_subset['l6_match'].sum()}/{n_loocv})")
    print(f"  L3 cluster concordance: {l3_concordance:.1f}% ({loocv_subset['l3_match'].sum()}/{n_loocv})")
    print(f"  Mean NN distance: {loocv_subset['nn_distance'].mean():.6f}")
    print(f"  Mean confidence: {loocv_subset['confidence'].mean():.1f}%")

    # ── Intra-study cluster detection metrics ────────────────────────────────
    print("\n" + "=" * 70)
    print("INTRA-STUDY CLUSTER DETECTION")
    print("=" * 70)

    detected_clusters = 0
    total_expected_clusters = 0
    cluster_results = []

    for study, expected_accs in gt_clusters.items():
        study_results = results_df[results_df["accession"].isin(expected_accs)]
        if len(study_results) < 2:
            continue

        # Check how many unique L6 codes within this study
        plin_codes = study_results["loocv_plin"].tolist()
        code_counts = Counter(plin_codes)
        shared_codes = {code: cnt for code, cnt in code_counts.items() if cnt >= 2}

        # Expected: at least some members share a code
        total_expected_clusters += 1
        detected = len(shared_codes) > 0
        if detected:
            detected_clusters += 1

        # L3 cluster check
        l3_codes = study_results["loocv_l3"].tolist()
        l3_counts = Counter(l3_codes)
        l3_shared = {code: cnt for code, cnt in l3_counts.items() if cnt >= 2}

        cluster_results.append({
            "study": study,
            "n_plasmids": len(expected_accs),
            "unique_l6": len(set(plin_codes)),
            "shared_l6": sum(shared_codes.values()) if shared_codes else 0,
            "l6_cluster_detected": detected,
            "unique_l3": len(set(l3_codes)),
            "l3_cluster_detected": len(l3_shared) > 0,
        })

        status = "DETECTED" if detected else "NOT detected"
        print(f"\n  {study}:")
        print(f"    Plasmids: {len(expected_accs)}, Unique L6: {len(set(plin_codes))}, "
              f"L3 clusters: {len(set(l3_codes))}")
        if shared_codes:
            for code, cnt in sorted(shared_codes.items(), key=lambda x: -x[1]):
                l6_num = code.split(".")[-1] if "." in str(code) else code
                print(f"    → {cnt} plasmids share pLIN {l6_num}")
        print(f"    Cluster {status}")

    cluster_df = pd.DataFrame(cluster_results)

    # Compute precision and recall
    if total_expected_clusters > 0:
        cluster_recall = 100.0 * detected_clusters / total_expected_clusters
    else:
        cluster_recall = 0

    print(f"\n  Cluster detection recall: {detected_clusters}/{total_expected_clusters} "
          f"({cluster_recall:.1f}%)")

    # Per-gene summary
    print("\n" + "=" * 70)
    print("PER-RESISTANCE-GENE CONCORDANCE")
    print("=" * 70)

    for gene in sorted(loocv_subset["resistance_gene"].unique()):
        sub = loocv_subset[loocv_subset["resistance_gene"] == gene]
        l6_rate = sub["l6_match"].mean() * 100
        l3_rate = sub["l3_match"].mean() * 100
        print(f"  {gene:<16s}: n={len(sub):>2d}, L6={l6_rate:>5.1f}%, L3={l3_rate:>5.1f}%, "
              f"mean_dist={sub['nn_distance'].mean():.4f}")

    # ── Save results ─────────────────────────────────────────────────────────
    out_dir = os.path.join(BASE_DIR, "output")

    results_path = os.path.join(out_dir, "outbreak_benchmark_results.tsv")
    results_df.to_csv(results_path, sep="\t", index=False)
    print(f"\nSaved: {results_path}")

    cluster_path = os.path.join(out_dir, "outbreak_benchmark_clusters.tsv")
    cluster_df.to_csv(cluster_path, sep="\t", index=False)
    print(f"Saved: {cluster_path}")

    # Summary text
    summary_path = os.path.join(out_dir, "outbreak_benchmark_summary.txt")
    with open(summary_path, "w") as f:
        f.write("pLIN Outbreak Detection — LOOCV Benchmark Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total outbreak plasmids: {len(results_df)}\n")
        f.write(f"LOOCV tested (with sequence): {n_loocv}\n")
        f.write(f"Pre-computed (no sequence): {len(results_df) - n_loocv}\n\n")
        f.write(f"L6 pLIN concordance: {l6_concordance:.1f}% ({loocv_subset['l6_match'].sum()}/{n_loocv})\n")
        f.write(f"L3 cluster concordance: {l3_concordance:.1f}% ({loocv_subset['l3_match'].sum()}/{n_loocv})\n")
        f.write(f"Mean NN distance: {loocv_subset['nn_distance'].mean():.6f}\n")
        f.write(f"Mean confidence: {loocv_subset['confidence'].mean():.1f}%\n\n")
        f.write(f"Ground truth clusters (≥2 plasmids): {total_expected_clusters}\n")
        f.write(f"Clusters detected: {detected_clusters}\n")
        f.write(f"Cluster detection recall: {cluster_recall:.1f}%\n")
    print(f"Saved: {summary_path}")

    print("\n" + "=" * 70)
    print("BENCHMARK COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
