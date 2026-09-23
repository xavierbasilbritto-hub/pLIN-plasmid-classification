#!/usr/bin/env python3
"""Compile combined outbreak validation results (original 17 + expanded 57 = 74 total)."""

import os
import pandas as pd

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Load original 17
df_orig = pd.read_csv(os.path.join(BASE_DIR, "output", "outbreak_validation_pLIN_results.tsv"), sep="\t")
df_orig["source"] = "original"
df_orig["country"] = df_orig["study"].apply(lambda s: s.split("_")[-1])

# Load expanded 57
df_exp = pd.read_csv(os.path.join(BASE_DIR, "output", "outbreak_validation_expanded_results.tsv"), sep="\t")
df_exp["source"] = "expanded"

# Standardize column names
df_orig_std = df_orig.rename(columns={"resistance_gene": "resistance_gene"})

# Combine
common_cols = ["accession", "study", "resistance_gene", "country", "length_bp",
               "predicted_inc", "confidence", "pLIN", "nn_plasmid", "nn_distance",
               "nn_inc_type", "nn_plin", "source"]

# Ensure both have the needed columns
for c in common_cols:
    if c not in df_orig_std.columns:
        df_orig_std[c] = ""
    if c not in df_exp.columns:
        df_exp[c] = ""

df_combined = pd.concat([df_orig_std[common_cols], df_exp[common_cols]], ignore_index=True)
df_combined.to_csv(os.path.join(BASE_DIR, "output", "outbreak_validation_combined_results.tsv"),
                   sep="\t", index=False)

# ── Statistics ────────────────────────────────────────────────────────────────
n = len(df_combined)
high_conf = df_combined[df_combined["confidence"] >= 60]
low_conf = df_combined[df_combined["confidence"] < 60]

print("=" * 70)
print("COMBINED OUTBREAK VALIDATION STATISTICS")
print("=" * 70)
print(f"  Total plasmids:           {n}")
print(f"  Original (4 studies):     {len(df_orig)}")
print(f"  Expanded (new):           {len(df_exp)}")
print(f"  Countries:                {df_combined['country'].nunique()}")
print(f"  Unique studies:           {df_combined['study'].nunique()}")
print(f"  Resistance genes:         {df_combined['resistance_gene'].nunique()}")
print(f"  Unique pLIN codes:        {df_combined['pLIN'].nunique()}")
print()
print(f"  High confidence (≥60%):   {len(high_conf)} ({100*len(high_conf)/n:.1f}%)")
print(f"  Low confidence (<60%):    {len(low_conf)} ({100*len(low_conf)/n:.1f}%)")
print()
print(f"  Mean NN distance:         {df_combined['nn_distance'].mean():.6f}")
print(f"  Median NN distance:       {df_combined['nn_distance'].median():.6f}")

# Count intra-study clusters
print("\n  Intra-study clustering (same pLIN within same study):")
cluster_count = 0
for study in sorted(df_combined["study"].unique()):
    sub = df_combined[df_combined["study"] == study]
    if len(sub) < 2:
        continue
    plin_counts = sub["pLIN"].value_counts()
    shared = plin_counts[plin_counts > 1]
    if not shared.empty:
        for code, cnt in shared.items():
            cluster_count += 1
            print(f"    {study}: {cnt} plasmids share pLIN {code}")

print(f"\n  Total intra-study clusters: {cluster_count}")

# Per-resistance-gene summary
print("\n  Per-gene summary:")
for gene in sorted(df_combined["resistance_gene"].unique()):
    sub = df_combined[df_combined["resistance_gene"] == gene]
    hc = sub[sub["confidence"] >= 60]
    print(f"    {gene:<16s}: n={len(sub):>2d}, high_conf={len(hc):>2d}/{len(sub):>2d} "
          f"({100*len(hc)/len(sub):.0f}%), mean_dist={sub['nn_distance'].mean():.4f}, "
          f"unique_pLIN={sub['pLIN'].nunique()}")

print("\n" + "=" * 70)
print("Done!")
