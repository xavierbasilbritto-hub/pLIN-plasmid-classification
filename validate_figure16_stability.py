#!/usr/bin/env python3
"""Compute real cluster stability and linkage comparison data for Figure 16.

Runs assess_cluster_stability() and compare_linkage_methods() on representative
subsets of training data per Inc group (to keep runtime feasible).
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "output", "figure16_validation")
os.makedirs(OUT_DIR, exist_ok=True)

CLASSIFIER_PATH = os.path.join(BASE_DIR, "data", "inc_classifier.npz")

PLIN_THRESHOLDS = {
    "L1": 0.150,
    "L2": 0.100,
    "L3": 0.050,
    "L4": 0.020,
    "L5": 0.010,
    "L6": 0.001,
}


def compute_cluster_stability_per_group():
    """Run bootstrap stability on each Inc group individually."""
    print("=" * 60)
    print("Panel A: Bootstrap cluster stability per group")
    print("=" * 60)

    data = np.load(CLASSIFIER_PATH, allow_pickle=True)
    X = data["X"].astype(np.float64)
    y = data["y"]
    group_names = list(data["group_names"])

    n_bootstrap = 50
    rng = np.random.default_rng(42)
    all_results = []

    for gi, gname in enumerate(group_names):
        mask = y == gi
        X_group = X[mask]
        n_group = len(X_group)

        if n_group < 4:
            print(f"  {gname}: n={n_group} (too small, skipping)")
            continue

        # Subsample large groups for feasibility
        if n_group > 200:
            idx_sub = rng.choice(n_group, 200, replace=False)
            X_sub = X_group[idx_sub]
            n_sub = 200
        else:
            X_sub = X_group
            n_sub = n_group

        # Original clustering
        dist_orig = pdist(X_sub, metric="cosine")
        Z_orig = linkage(dist_orig, method="single")

        # Get clusters at L3 (the most informative level for stability)
        for level_name, thresh in [("L3", 0.050), ("L5", 0.010)]:
            clusters_orig = fcluster(Z_orig, t=thresh, criterion="distance")
            n_clusters = len(set(clusters_orig))

            if n_clusters < 2:
                all_results.append({
                    "inc_group": gname,
                    "level": level_name,
                    "n_samples": n_sub,
                    "n_clusters": n_clusters,
                    "mean_stability": 100.0,
                    "min_stability": 100.0,
                    "n_stable": n_clusters,
                    "n_moderate": 0,
                    "n_unstable": 0,
                })
                continue

            # Bootstrap: track co-clustering
            co_cluster = np.zeros((n_sub, n_sub))
            for b in range(n_bootstrap):
                boot_idx = rng.choice(n_sub, n_sub, replace=True)
                boot_vectors = X_sub[boot_idx]
                boot_dist = pdist(boot_vectors, metric="cosine")
                Z_boot = linkage(boot_dist, method="single")
                boot_clusters = fcluster(Z_boot, t=thresh, criterion="distance")

                for i in range(n_sub):
                    for j in range(i + 1, n_sub):
                        if boot_clusters[i] == boot_clusters[j]:
                            co_cluster[boot_idx[i], boot_idx[j]] += 1
                            co_cluster[boot_idx[j], boot_idx[i]] += 1

            # Compute per-cluster stability
            cluster_stabilities = []
            for cl_id in set(clusters_orig):
                members = [i for i, c in enumerate(clusters_orig) if c == cl_id]
                if len(members) < 2:
                    cluster_stabilities.append(100.0)
                    continue
                min_support = 100.0
                for i in range(len(members)):
                    for j in range(i + 1, len(members)):
                        support = 100.0 * co_cluster[members[i], members[j]] / n_bootstrap
                        min_support = min(min_support, support)
                cluster_stabilities.append(min_support)

            stab_arr = np.array(cluster_stabilities)
            n_stable = int(np.sum(stab_arr >= 80))
            n_moderate = int(np.sum((stab_arr >= 50) & (stab_arr < 80)))
            n_unstable = int(np.sum(stab_arr < 50))

            all_results.append({
                "inc_group": gname,
                "level": level_name,
                "n_samples": n_sub,
                "n_clusters": n_clusters,
                "mean_stability": round(float(np.mean(stab_arr)), 1),
                "min_stability": round(float(np.min(stab_arr)), 1),
                "n_stable": n_stable,
                "n_moderate": n_moderate,
                "n_unstable": n_unstable,
            })

        print(f"  {gname}: n={n_sub}, done")

    df = pd.DataFrame(all_results)
    out_path = os.path.join(OUT_DIR, "cluster_stability.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\n  Saved: {out_path}")
    print(f"  Total entries: {len(df)}")

    return df


def compute_linkage_comparison():
    """Run linkage comparison on representative groups."""
    from sklearn.metrics import adjusted_rand_score

    print("\n" + "=" * 60)
    print("Panel B: Linkage method comparison")
    print("=" * 60)

    data = np.load(CLASSIFIER_PATH, allow_pickle=True)
    X = data["X"].astype(np.float64)
    y = data["y"]
    group_names = list(data["group_names"])

    rng = np.random.default_rng(42)
    all_results = []

    # Run on each group
    for gi, gname in enumerate(group_names):
        mask = y == gi
        X_group = X[mask]
        n_group = len(X_group)

        if n_group < 10:
            continue

        if n_group > 200:
            idx_sub = rng.choice(n_group, 200, replace=False)
            X_sub = X_group[idx_sub]
        else:
            X_sub = X_group

        dist_condensed = pdist(X_sub, metric="cosine")

        methods = ["single", "complete", "average"]
        clusters_by_method = {}
        for method in methods:
            Z = linkage(dist_condensed, method=method)
            clusters_by_method[method] = {}
            for level, thresh in PLIN_THRESHOLDS.items():
                clusters_by_method[method][level] = fcluster(Z, t=thresh, criterion="distance")

        # ARI at L3 level (most discriminating)
        for level in ["L3", "L5"]:
            s = clusters_by_method["single"][level]
            c = clusters_by_method["complete"][level]
            a = clusters_by_method["average"][level]

            all_results.append({
                "inc_group": gname,
                "level": level,
                "n_samples": len(X_sub),
                "single_vs_complete": round(adjusted_rand_score(s, c), 3),
                "single_vs_average": round(adjusted_rand_score(s, a), 3),
                "complete_vs_average": round(adjusted_rand_score(c, a), 3),
                "n_clusters_single": len(set(s)),
                "n_clusters_complete": len(set(c)),
                "n_clusters_average": len(set(a)),
            })

        print(f"  {gname}: n={len(X_sub)}, done")

    df = pd.DataFrame(all_results)
    out_path = os.path.join(OUT_DIR, "linkage_comparison.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\n  Saved: {out_path}")

    # Summary: mean ARI across all groups at L3
    l3 = df[df["level"] == "L3"]
    print(f"\n  L3 mean ARI (single vs complete): {l3['single_vs_complete'].mean():.3f}")
    print(f"  L3 mean ARI (single vs average):  {l3['single_vs_average'].mean():.3f}")
    print(f"  L3 mean ARI (complete vs average): {l3['complete_vs_average'].mean():.3f}")
    print(f"  L3 mean clusters (single):   {l3['n_clusters_single'].mean():.1f}")
    print(f"  L3 mean clusters (complete):  {l3['n_clusters_complete'].mean():.1f}")
    print(f"  L3 mean clusters (average):   {l3['n_clusters_average'].mean():.1f}")

    return df


if __name__ == "__main__":
    stab_df = compute_cluster_stability_per_group()
    link_df = compute_linkage_comparison()

    # Save combined summary
    l3_stab = stab_df[stab_df["level"] == "L3"]
    l3_link = link_df[link_df["level"] == "L3"]

    # Collect all stability scores for the histogram
    all_stab_scores = []
    for _, row in l3_stab.iterrows():
        # We stored aggregates per group; for the histogram we need individual scores
        all_stab_scores.extend([row["mean_stability"]] * row["n_clusters"])

    summary = {
        "stability": {
            "n_groups_assessed": len(l3_stab),
            "total_clusters_L3": int(l3_stab["n_clusters"].sum()),
            "total_stable_L3": int(l3_stab["n_stable"].sum()),
            "total_moderate_L3": int(l3_stab["n_moderate"].sum()),
            "total_unstable_L3": int(l3_stab["n_unstable"].sum()),
            "overall_mean_stability": round(float(l3_stab["mean_stability"].mean()), 1),
        },
        "linkage": {
            "n_groups_assessed": len(l3_link),
            "L3_mean_ari_single_vs_complete": round(float(l3_link["single_vs_complete"].mean()), 3),
            "L3_mean_ari_single_vs_average": round(float(l3_link["single_vs_average"].mean()), 3),
            "L3_mean_ari_complete_vs_average": round(float(l3_link["complete_vs_average"].mean()), 3),
            "L3_mean_n_clusters_single": round(float(l3_link["n_clusters_single"].mean()), 1),
            "L3_mean_n_clusters_complete": round(float(l3_link["n_clusters_complete"].mean()), 1),
            "L3_mean_n_clusters_average": round(float(l3_link["n_clusters_average"].mean()), 1),
        }
    }

    summary_path = os.path.join(OUT_DIR, "figure16_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Summary: {summary_path}")
    print("\nDone!")
