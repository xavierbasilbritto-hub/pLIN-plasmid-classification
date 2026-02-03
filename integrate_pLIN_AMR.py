#!/usr/bin/env python3
"""
pLIN + AMRFinderPlus Integration Script
Cross-references pLIN hierarchical codes with AMRFinderPlus gene detections.
Produces a combined surveillance table and summary statistics.
"""

import os
import pandas as pd
import numpy as np
from collections import Counter, defaultdict

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PLIN_FILE = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")
AMR_FILE = os.path.join(BASE_DIR, "output", "amrfinder", "amrfinder_all_plasmids.tsv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "integrated")


def load_data():
    """Load pLIN assignments and AMRFinderPlus results."""
    plin = pd.read_csv(PLIN_FILE, sep="\t")
    print(f"pLIN assignments loaded: {len(plin)} plasmids, {plin['pLIN'].nunique()} unique codes")

    amr = pd.read_csv(AMR_FILE, sep="\t")
    print(f"AMRFinderPlus hits loaded: {len(amr)} gene detections")
    print(f"  Plasmids with hits: {amr['source_file'].nunique()}")
    print(f"  Element types: {amr['Type'].value_counts().to_dict()}")
    return plin, amr


def build_per_plasmid_amr_summary(amr):
    """Summarise AMR/virulence genes per plasmid."""
    records = []
    for source, grp in amr.groupby("source_file"):
        # Separate AMR from virulence/stress
        amr_genes = grp[grp["Type"] == "AMR"]
        vir_genes = grp[grp["Type"] == "VIRULENCE"]
        stress_genes = grp[grp["Type"] == "STRESS"]

        amr_symbols = sorted(amr_genes["Element symbol"].unique())
        vir_symbols = sorted(vir_genes["Element symbol"].unique())
        stress_symbols = sorted(stress_genes["Element symbol"].unique())

        amr_classes = sorted(amr_genes["Class"].dropna().unique())
        amr_subclasses = sorted(amr_genes["Subclass"].dropna().unique())

        records.append({
            "source_file": source,
            "n_amr_genes": len(amr_genes),
            "n_vir_genes": len(vir_genes),
            "n_stress_genes": len(stress_genes),
            "n_total_hits": len(grp),
            "amr_genes": "; ".join(amr_symbols) if amr_symbols else "none",
            "vir_genes": "; ".join(vir_symbols) if vir_symbols else "none",
            "stress_genes": "; ".join(stress_symbols) if stress_symbols else "none",
            "amr_classes": "; ".join(amr_classes) if amr_classes else "none",
            "amr_subclasses": "; ".join(amr_subclasses) if amr_subclasses else "none",
            "inc_type": grp["inc_type"].iloc[0],
        })

    return pd.DataFrame(records)


def merge_plin_amr(plin, amr_summary):
    """Merge pLIN codes with AMR summary per plasmid."""
    # Match on plasmid_id (pLIN) ↔ source_file (AMR)
    # source_file is like "RefSeq_NC_002120.1", plasmid_id is "RefSeq_NC_002120.1"
    merged = plin.merge(amr_summary, left_on="plasmid_id", right_on="source_file", how="left")

    # Fill missing (plasmids with no AMR hits)
    fill_cols = ["n_amr_genes", "n_vir_genes", "n_stress_genes", "n_total_hits"]
    for col in fill_cols:
        merged[col] = merged[col].fillna(0).astype(int)
    for col in ["amr_genes", "vir_genes", "stress_genes", "amr_classes", "amr_subclasses"]:
        merged[col] = merged[col].fillna("none")

    return merged


def analyse_plin_amr(merged):
    """Produce statistics on AMR gene distribution across pLIN lineages."""
    print("\n" + "=" * 80)
    print("pLIN + AMRFinderPlus INTEGRATED ANALYSIS")
    print("=" * 80)

    # ── Overview ──
    total = len(merged)
    with_amr = (merged["n_amr_genes"] > 0).sum()
    with_vir = (merged["n_vir_genes"] > 0).sum()
    with_stress = (merged["n_stress_genes"] > 0).sum()
    with_any = (merged["n_total_hits"] > 0).sum()

    print(f"\n1. OVERALL AMR/VIRULENCE PREVALENCE")
    print(f"   Total plasmids:              {total}")
    print(f"   Plasmids with AMR genes:     {with_amr} ({with_amr/total*100:.1f}%)")
    print(f"   Plasmids with virulence:     {with_vir} ({with_vir/total*100:.1f}%)")
    print(f"   Plasmids with stress genes:  {with_stress} ({with_stress/total*100:.1f}%)")
    print(f"   Plasmids with any hit:       {with_any} ({with_any/total*100:.1f}%)")

    # ── Per Inc type ──
    print(f"\n2. AMR PREVALENCE BY INC TYPE")
    print(f"   {'Inc Type':<10s} {'Total':>6s} {'AMR+':>6s} {'%AMR':>7s} {'VIR+':>6s} {'%VIR':>7s} {'Mean AMR genes':>15s}")
    print(f"   {'-'*60}")
    for inc in ["IncFII", "IncN", "IncX1"]:
        sub = merged[merged["inc_type_x"] == inc]
        n = len(sub)
        n_amr = (sub["n_amr_genes"] > 0).sum()
        n_vir = (sub["n_vir_genes"] > 0).sum()
        mean_amr = sub["n_amr_genes"].mean()
        print(f"   {inc:<10s} {n:>6d} {n_amr:>6d} {n_amr/n*100:>6.1f}% {n_vir:>6d} {n_vir/n*100:>6.1f}% {mean_amr:>14.2f}")

    # ── Most common AMR genes ──
    print(f"\n3. MOST COMMON AMR GENES (Top 20)")
    amr_only = merged[merged["n_amr_genes"] > 0]
    all_amr_genes = []
    for genes_str in amr_only["amr_genes"]:
        if genes_str != "none":
            all_amr_genes.extend(genes_str.split("; "))
    gene_counts = Counter(all_amr_genes)
    print(f"   {'Gene':<25s} {'Count':>6s} {'% of AMR+ plasmids':>20s}")
    print(f"   {'-'*55}")
    n_amr_plasmids = len(amr_only)
    for gene, count in gene_counts.most_common(20):
        print(f"   {gene:<25s} {count:>6d} {count/n_amr_plasmids*100:>19.1f}%")

    # ── Most common AMR classes ──
    print(f"\n4. AMR DRUG CLASSES")
    all_classes = []
    for c in merged["amr_classes"]:
        if c != "none":
            all_classes.extend(c.split("; "))
    class_counts = Counter(all_classes)
    print(f"   {'Drug Class':<40s} {'Plasmids':>8s}")
    print(f"   {'-'*50}")
    for cls, count in class_counts.most_common(15):
        print(f"   {cls:<40s} {count:>8d}")

    # ── pLIN lineages as AMR vehicles ──
    print(f"\n5. TOP pLIN LINEAGES CARRYING AMR GENES")
    plin_amr = merged[merged["n_amr_genes"] > 0].groupby("pLIN").agg(
        n_plasmids=("plasmid_id", "count"),
        mean_amr=("n_amr_genes", "mean"),
        total_amr=("n_amr_genes", "sum"),
        inc_types=("inc_type_x", lambda x: ", ".join(sorted(x.unique()))),
    ).sort_values("n_plasmids", ascending=False)

    print(f"   {'pLIN Code':<30s} {'n':>5s} {'Mean AMR':>9s} {'Total AMR':>10s} {'Inc Type(s)'}")
    print(f"   {'-'*80}")
    for code, row in plin_amr.head(20).iterrows():
        print(f"   {code:<30s} {row['n_plasmids']:>5d} {row['mean_amr']:>8.1f} {row['total_amr']:>10.0f} {row['inc_types']}")

    # ── AMR gene profiles per pLIN code ──
    print(f"\n6. AMR GENE PROFILES OF TOP LINEAGES")
    for code in plin_amr.head(10).index:
        sub = merged[(merged["pLIN"] == code) & (merged["n_amr_genes"] > 0)]
        genes = []
        for g in sub["amr_genes"]:
            if g != "none":
                genes.extend(g.split("; "))
        gene_freq = Counter(genes)
        n = len(sub)
        top_genes = [f"{g} ({c}/{n})" for g, c in gene_freq.most_common(5)]
        print(f"   pLIN {code}: n={n}, top genes: {', '.join(top_genes)}")

    # ── Clinically critical genes ──
    print(f"\n7. CLINICALLY CRITICAL RESISTANCE GENES")
    critical_patterns = {
        "Carbapenemases": ["blaKPC", "blaNDM", "blaOXA-48", "blaVIM", "blaIMP"],
        "ESBLs": ["blaCTX-M", "blaSHV", "blaTEM"],
        "Colistin resistance": ["mcr-"],
        "Vancomycin resistance": ["van"],
        "PMQR (fluoroquinolone)": ["qnr", "aac(6')-Ib-cr", "oqxA", "oqxB"],
    }

    all_genes_flat = []
    for genes_str in merged["amr_genes"]:
        if genes_str != "none":
            all_genes_flat.extend(genes_str.split("; "))

    for category, patterns in critical_patterns.items():
        matching = [g for g in all_genes_flat if any(p in g for p in patterns)]
        if matching:
            gene_counts_crit = Counter(matching)
            total_hits = sum(gene_counts_crit.values())
            top = ", ".join(f"{g}({c})" for g, c in gene_counts_crit.most_common(5))
            print(f"   {category}: {total_hits} detections — {top}")
        else:
            print(f"   {category}: none detected")

    # ── Virulence gene distribution ──
    print(f"\n8. VIRULENCE GENE PREVALENCE BY INC TYPE")
    for inc in ["IncFII", "IncN", "IncX1"]:
        sub = merged[(merged["inc_type_x"] == inc) & (merged["n_vir_genes"] > 0)]
        vir_genes_all = []
        for g in sub["vir_genes"]:
            if g != "none":
                vir_genes_all.extend(g.split("; "))
        vir_counts = Counter(vir_genes_all)
        top5 = ", ".join(f"{g}({c})" for g, c in vir_counts.most_common(5))
        print(f"   {inc}: {len(sub)} plasmids with virulence genes — top: {top5}")

    return plin_amr


def main():
    print("=" * 80)
    print("pLIN + AMRFinderPlus Integration")
    print("=" * 80)

    print("\nLoading data ...")
    plin, amr = load_data()

    print("\nBuilding per-plasmid AMR summary ...")
    amr_summary = build_per_plasmid_amr_summary(amr)
    print(f"  {len(amr_summary)} plasmids with AMRFinderPlus detections")

    print("\nMerging pLIN codes with AMR data ...")
    merged = merge_plin_amr(plin, amr_summary)

    # Save combined table
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    out_file = os.path.join(OUTPUT_DIR, "pLIN_AMR_integrated.tsv")
    cols_to_save = [
        "plasmid_id", "inc_type_x", "length_bp", "pLIN",
        "bin_A", "bin_B", "bin_C", "bin_D", "bin_E", "bin_F",
        "n_amr_genes", "n_vir_genes", "n_stress_genes", "n_total_hits",
        "amr_genes", "vir_genes", "stress_genes", "amr_classes", "amr_subclasses",
    ]
    merged[cols_to_save].to_csv(out_file, sep="\t", index=False)
    print(f"  Integrated table saved: {out_file}")

    # Analyse
    plin_amr = analyse_plin_amr(merged)

    # Save pLIN lineage AMR summary
    lineage_file = os.path.join(OUTPUT_DIR, "pLIN_lineage_AMR_summary.tsv")
    plin_amr.to_csv(lineage_file, sep="\t")
    print(f"\n  Lineage AMR summary saved: {lineage_file}")

    print("\n" + "=" * 80)
    print("Integration complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
