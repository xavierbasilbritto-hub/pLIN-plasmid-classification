#!/usr/bin/env python3
"""
pLIN Test — Integrate AMR + Generate Cladograms with AMR annotations
1. Integrates pLIN assignments with AMRFinderPlus results for test plasmids
2. Generates cladograms annotated with AMR gene presence/absence
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap
from itertools import product as iter_product
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from Bio import SeqIO

# ── Configuration ──────────────────────────────────────────────────────────────

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.join(BASE_DIR, "test_plasmids")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "test")

PLIN_FILE = os.path.join(OUTPUT_DIR, "test_pLIN_assignments.tsv")
AMR_FILE = os.path.join(OUTPUT_DIR, "amrfinder", "amrfinder_test_all.tsv")

PLIN_THRESHOLDS = {
    "A": 0.150, "B": 0.100, "C": 0.050,
    "D": 0.020, "E": 0.010, "F": 0.001,
}

THRESHOLD_COLORS = {
    "A": "#E53935", "B": "#FB8C00", "C": "#FDD835",
    "D": "#43A047", "E": "#1E88E5", "F": "#8E24AA",
}

PLIN_LEVEL_NAMES = {
    "A": "Family", "B": "Subfamily", "C": "Cluster",
    "D": "Subcluster", "E": "Clone", "F": "Strain",
}

STRAIN_COLORS = [
    "#2196F3", "#FF9800", "#4CAF50", "#E91E63", "#9C27B0",
    "#00BCD4", "#FF5722", "#795548", "#607D8B", "#CDDC39",
    "#F44336", "#3F51B5", "#009688", "#FFC107", "#8BC34A",
]

TYPE_COLORS = {
    "AMR": "#E53935",
    "STRESS": "#FB8C00",
    "VIRULENCE": "#8E24AA",
}


# ── Data loading ─────────────────────────────────────────────────────────────

def load_sequences():
    records = []
    for entry in sorted(os.listdir(TEST_DIR)):
        subdir = os.path.join(TEST_DIR, entry)
        if not os.path.isdir(subdir):
            continue
        fasta_files = sorted(glob.glob(os.path.join(subdir, "*.fasta")))
        if not fasta_files:
            continue
        for fasta_file in fasta_files:
            for rec in SeqIO.parse(fasta_file, "fasta"):
                records.append({
                    "plasmid_id": rec.id,
                    "inc_type": entry,
                    "sequence": str(rec.seq),
                    "length": len(rec.seq),
                    "source_file": os.path.basename(fasta_file),
                })
    return records


def compute_kmer_vectors(records, k=4):
    bases = "ACGT"
    all_kmers = ["".join(p) for p in iter_product(bases, repeat=k)]
    vectors = np.zeros((len(records), len(all_kmers)), dtype=np.float64)
    for idx, rec in enumerate(records):
        seq = rec["sequence"].upper()
        total = max(len(seq) - k + 1, 1)
        for ki, kmer in enumerate(all_kmers):
            vectors[idx, ki] = seq.count(kmer) / total
    return vectors


def get_plin_codes(records, Z):
    codes = {}
    for bname, thresh in PLIN_THRESHOLDS.items():
        codes[bname] = fcluster(Z, t=thresh, criterion="distance")
    plin = []
    for i in range(len(records)):
        parts = [str(codes[b][i]) for b in PLIN_THRESHOLDS]
        plin.append(".".join(parts))
    return plin, codes


def get_strain_cmap(strain_clusters):
    unique_strains = sorted(set(strain_clusters))
    return {s: STRAIN_COLORS[i % len(STRAIN_COLORS)] for i, s in enumerate(unique_strains)}


# ── Part 1: Integrate pLIN + AMR ────────────────────────────────────────────

def integrate_plin_amr():
    print("=" * 70)
    print("PART 1: Integrating pLIN + AMRFinderPlus for test plasmids")
    print("=" * 70)

    plin_df = pd.read_csv(PLIN_FILE, sep="\t")
    amr_df = pd.read_csv(AMR_FILE, sep="\t")

    print(f"\n  pLIN assignments: {len(plin_df)} plasmids")
    print(f"  AMR detections:   {len(amr_df)} total hits")

    # Map source_file to match pLIN source_file (add .fasta)
    amr_df["source_file_fasta"] = amr_df["source_file"] + ".fasta"

    # Summary per plasmid
    plasmid_amr = []
    for _, row in plin_df.iterrows():
        sf = row["source_file"] if "source_file" in plin_df.columns else row["plasmid_id"]
        # Match by source file name without .fasta
        sf_base = sf.replace(".fasta", "") if isinstance(sf, str) else sf
        hits = amr_df[amr_df["source_file"] == sf_base]

        amr_hits = hits[hits["Type"] == "AMR"]
        stress_hits = hits[hits["Type"] == "STRESS"]
        vir_hits = hits[hits["Type"] == "VIRULENCE"]

        amr_genes = sorted(amr_hits["Element symbol"].unique()) if len(amr_hits) > 0 else []
        stress_genes = sorted(stress_hits["Element symbol"].unique()) if len(stress_hits) > 0 else []

        plasmid_amr.append({
            "plasmid_id": row["plasmid_id"],
            "inc_type": row["inc_type"],
            "length_bp": row["length_bp"],
            "pLIN": row["pLIN"],
            "total_hits": len(hits),
            "AMR_count": len(amr_hits),
            "STRESS_count": len(stress_hits),
            "VIRULENCE_count": len(vir_hits),
            "AMR_genes": "; ".join(amr_genes),
            "STRESS_genes": "; ".join(stress_genes),
            "AMR_classes": "; ".join(sorted(amr_hits["Class"].unique())) if len(amr_hits) > 0 else "",
        })

    integrated = pd.DataFrame(plasmid_amr)
    out_path = os.path.join(OUTPUT_DIR, "test_pLIN_AMR_integrated.tsv")
    integrated.to_csv(out_path, sep="\t", index=False)
    print(f"\n  Integrated table saved: {out_path}")

    # Summary statistics
    print(f"\n  --- Summary ---")
    print(f"  Plasmids with AMR genes: {(integrated['AMR_count'] > 0).sum()}/{len(integrated)}")
    print(f"  Plasmids with STRESS genes: {(integrated['STRESS_count'] > 0).sum()}/{len(integrated)}")
    print(f"  Total AMR detections: {integrated['AMR_count'].sum()}")
    print(f"  Total STRESS detections: {integrated['STRESS_count'].sum()}")

    print(f"\n  AMR gene distribution:")
    amr_only = amr_df[amr_df["Type"] == "AMR"]
    for gene, count in amr_only["Element symbol"].value_counts().items():
        n_plasmids = amr_only[amr_only["Element symbol"] == gene]["source_file"].nunique()
        print(f"    {gene:<15s}  {count:>3d} detections in {n_plasmids:>2d} plasmids")

    print(f"\n  AMR drug classes:")
    for cls, count in amr_only["Class"].value_counts().items():
        print(f"    {cls:<20s}  {count:>3d}")

    return integrated, amr_df


# ── Part 2: Cladogram with AMR annotations ──────────────────────────────────

def plot_cladogram_amr(Z, labels, plin_codes, strain_clusters, integrated, amr_df, records):
    """Rectangular cladogram with AMR gene presence/absence heatmap."""
    strain_cmap = get_strain_cmap(strain_clusters)
    unique_strains = sorted(set(strain_clusters))

    # Get key AMR/stress genes for heatmap
    all_genes = sorted(amr_df["Element symbol"].unique())
    # Group: AMR genes first, then stress
    amr_genes = sorted(amr_df[amr_df["Type"] == "AMR"]["Element symbol"].unique())
    stress_genes = sorted(amr_df[amr_df["Type"] == "STRESS"]["Element symbol"].unique())
    gene_list = amr_genes + stress_genes
    gene_types = ["AMR"] * len(amr_genes) + ["STRESS"] * len(stress_genes)

    n_genes = len(gene_list)
    n_plasmids = len(labels)

    # Build presence/absence matrix
    label_to_source = {labels[i]: records[i]["source_file"].replace(".fasta", "") for i in range(n_plasmids)}

    fig = plt.figure(figsize=(20, 12))

    # Layout
    ax_dendro = fig.add_axes([0.01, 0.10, 0.18, 0.78])
    ax_plin = fig.add_axes([0.20, 0.10, 0.10, 0.78])
    ax_heat = fig.add_axes([0.32, 0.10, 0.45, 0.78])
    ax_legend = fig.add_axes([0.80, 0.10, 0.18, 0.78])

    # -- Dendrogram --
    with plt.rc_context({"lines.linewidth": 2.0}):
        ddata = dendrogram(
            Z, labels=labels, orientation="right",
            leaf_font_size=9, ax=ax_dendro,
            color_threshold=0, above_threshold_color="#444444",
        )

    max_merge = Z[:, 2].max()
    ax_dendro.set_xlim(0, max_merge * 1.15)

    leaf_colors = {labels[i]: strain_cmap[sc] for i, sc in enumerate(strain_clusters)}
    for lbl in ax_dendro.get_yticklabels():
        txt = lbl.get_text()
        if txt in leaf_colors:
            lbl.set_color(leaf_colors[txt])
            lbl.set_fontweight("bold")

    for bname, thresh in PLIN_THRESHOLDS.items():
        if thresh <= max_merge * 1.15:
            ax_dendro.axvline(x=thresh, color=THRESHOLD_COLORS[bname],
                              linestyle="--", linewidth=1.0, alpha=0.6)

    ax_dendro.set_xlabel("Cosine Distance", fontsize=9)
    ax_dendro.set_title("Dendrogram", fontsize=10, fontweight="bold")

    leaf_order = ddata["leaves"]

    # -- pLIN codes panel --
    ax_plin.set_ylim(ax_dendro.get_ylim())
    ax_plin.set_xlim(0, 1)
    ax_plin.axis("off")
    ax_plin.set_title("pLIN", fontsize=10, fontweight="bold")

    y_positions = [5 + i * 10 for i in range(len(leaf_order))]
    label_to_plin = dict(zip(labels, plin_codes))
    label_to_strain = dict(zip(labels, strain_clusters))

    for i, leaf_idx in enumerate(leaf_order):
        lbl = labels[leaf_idx]
        plin = label_to_plin[lbl]
        sc = label_to_strain[lbl]
        color = strain_cmap[sc]
        ax_plin.text(0.05, y_positions[i], plin, fontsize=8, fontfamily="monospace",
                     fontweight="bold", color=color, va="center")

    # -- AMR heatmap --
    heat_data = np.zeros((n_plasmids, n_genes))
    for i, leaf_idx in enumerate(leaf_order):
        source = label_to_source[labels[leaf_idx]]
        hits = amr_df[amr_df["source_file"] == source]
        for j, gene in enumerate(gene_list):
            if gene in hits["Element symbol"].values:
                heat_data[i, j] = 2 if gene_types[j] == "AMR" else 1

    # Custom colormap: white=absent, orange=stress, red=AMR
    cmap = ListedColormap(["#F5F5F5", "#FFE0B2", "#EF5350"])

    ax_heat.imshow(heat_data, aspect="auto", cmap=cmap, interpolation="nearest",
                   vmin=0, vmax=2)

    # Gene labels on top
    ax_heat.set_xticks(range(n_genes))
    ax_heat.set_xticklabels(gene_list, fontsize=7, rotation=65, ha="left",
                             rotation_mode="anchor")
    # Color gene labels by type
    for j, tick in enumerate(ax_heat.get_xticklabels()):
        tick.set_color(TYPE_COLORS.get(gene_types[j], "#333333"))
        tick.set_fontweight("bold")

    ax_heat.set_yticks(range(n_plasmids))
    ax_heat.set_yticklabels([labels[i] for i in leaf_order], fontsize=8)

    # Mark present cells
    for i in range(heat_data.shape[0]):
        for j in range(heat_data.shape[1]):
            if heat_data[i, j] > 0:
                ax_heat.text(j, i, "\u2713", ha="center", va="center",
                             fontsize=7, fontweight="bold",
                             color="white" if heat_data[i, j] == 2 else "#E65100")

    # Grid lines
    ax_heat.set_xticks(np.arange(-0.5, n_genes), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, n_plasmids), minor=True)
    ax_heat.grid(which="minor", color="#E0E0E0", linewidth=0.5)
    ax_heat.tick_params(which="minor", size=0)

    # Separator line between AMR and STRESS genes
    if amr_genes and stress_genes:
        sep_x = len(amr_genes) - 0.5
        ax_heat.axvline(x=sep_x, color="black", linewidth=2, linestyle="-")
        # Labels
        ax_heat.text(len(amr_genes) / 2 - 0.5, -1.5, "AMR genes",
                     ha="center", fontsize=9, fontweight="bold", color=TYPE_COLORS["AMR"])
        ax_heat.text(len(amr_genes) + len(stress_genes) / 2 - 0.5, -1.5, "Stress genes",
                     ha="center", fontsize=9, fontweight="bold", color=TYPE_COLORS["STRESS"])

    ax_heat.set_title("Gene Presence / Absence", fontsize=10, fontweight="bold")

    # -- Legend panel --
    ax_legend.axis("off")

    # Strain legend
    y = 0.95
    ax_legend.text(0.05, y, "pLIN Strain Clusters (F)", fontsize=9,
                   fontweight="bold", transform=ax_legend.transAxes)
    y -= 0.04
    for sc in unique_strains:
        n = sum(1 for s in strain_clusters if s == sc)
        ax_legend.add_patch(mpatches.FancyBboxPatch(
            (0.05, y - 0.01), 0.08, 0.025, transform=ax_legend.transAxes,
            boxstyle="round,pad=0.003", facecolor=strain_cmap[sc], edgecolor="none"))
        ax_legend.text(0.16, y + 0.003, f"Strain {sc} (n={n})", fontsize=8,
                       transform=ax_legend.transAxes, va="center")
        y -= 0.035

    # Type legend
    y -= 0.03
    ax_legend.text(0.05, y, "Detection Type", fontsize=9,
                   fontweight="bold", transform=ax_legend.transAxes)
    y -= 0.04
    for typ, color in TYPE_COLORS.items():
        if typ in amr_df["Type"].values:
            ax_legend.add_patch(mpatches.FancyBboxPatch(
                (0.05, y - 0.01), 0.08, 0.025, transform=ax_legend.transAxes,
                boxstyle="round,pad=0.003", facecolor=color if typ == "AMR" else "#FFE0B2",
                edgecolor="none"))
            ax_legend.text(0.16, y + 0.003, typ, fontsize=8,
                           transform=ax_legend.transAxes, va="center",
                           color=color, fontweight="bold")
            y -= 0.035

    # AMR summary
    y -= 0.03
    ax_legend.text(0.05, y, "AMR Summary", fontsize=9,
                   fontweight="bold", transform=ax_legend.transAxes)
    y -= 0.035

    amr_only = amr_df[amr_df["Type"] == "AMR"]
    for gene in amr_genes:
        n = amr_only[amr_only["Element symbol"] == gene]["source_file"].nunique()
        cls = amr_only[amr_only["Element symbol"] == gene]["Class"].iloc[0]
        ax_legend.text(0.05, y, f"\u2022 {gene} ({cls}) — {n} plasmids",
                       fontsize=7, transform=ax_legend.transAxes, va="center",
                       color=TYPE_COLORS["AMR"])
        y -= 0.03

    # Threshold info
    y -= 0.03
    ax_legend.text(0.05, y, "pLIN Thresholds", fontsize=9,
                   fontweight="bold", transform=ax_legend.transAxes)
    y -= 0.035
    for bname, thresh in PLIN_THRESHOLDS.items():
        ax_legend.text(0.05, y,
                       f"{bname} ({PLIN_LEVEL_NAMES[bname]}): d\u2264{thresh}",
                       fontsize=7, transform=ax_legend.transAxes, va="center",
                       color=THRESHOLD_COLORS[bname], fontweight="bold")
        y -= 0.025

    fig.suptitle(
        "pLIN Cladogram with AMR Gene Profile \u2014 IncX Test Plasmids (n=22)\n"
        "AMRFinderPlus v4.2.5 | 104 detections (34 AMR + 70 Stress)",
        fontsize=14, fontweight="bold", y=0.98,
    )

    return fig


def plot_amr_summary(integrated, amr_df):
    """Summary figure: AMR burden per plasmid and gene frequency."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel A: AMR + Stress gene count per plasmid
    ax = axes[0, 0]
    df = integrated.sort_values("total_hits", ascending=True)
    y = range(len(df))
    bars_amr = ax.barh(y, df["AMR_count"].values, color=TYPE_COLORS["AMR"],
                       label="AMR", height=0.7)
    bars_stress = ax.barh(y, df["STRESS_count"].values, left=df["AMR_count"].values,
                          color=TYPE_COLORS["STRESS"], label="Stress", height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{row['pLIN']}" for _, row in df.iterrows()], fontsize=7)
    ax.set_xlabel("Number of Detections", fontsize=10)
    ax.set_title("A. Gene Detections per Plasmid (by pLIN)", fontsize=11, fontweight="bold")
    ax.legend(fontsize=9)

    # Panel B: AMR gene frequency across plasmids
    ax = axes[0, 1]
    amr_only = amr_df[amr_df["Type"] == "AMR"]
    gene_freq = amr_only.groupby("Element symbol")["source_file"].nunique().sort_values(ascending=True)
    colors = [TYPE_COLORS["AMR"]] * len(gene_freq)
    ax.barh(range(len(gene_freq)), gene_freq.values, color=colors, height=0.7)
    ax.set_yticks(range(len(gene_freq)))
    ax.set_yticklabels(gene_freq.index, fontsize=9)
    ax.set_xlabel("Number of Plasmids", fontsize=10)
    ax.set_title("B. AMR Gene Prevalence", fontsize=11, fontweight="bold")

    # Panel C: Drug class distribution
    ax = axes[1, 0]
    class_counts = amr_only["Class"].value_counts()
    wedges, texts, autotexts = ax.pie(
        class_counts.values, labels=class_counts.index,
        autopct="%1.0f%%", startangle=90,
        colors=["#EF5350", "#42A5F5", "#66BB6A", "#FFA726"],
    )
    for t in autotexts:
        t.set_fontsize(10)
        t.set_fontweight("bold")
    ax.set_title("C. AMR Drug Classes", fontsize=11, fontweight="bold")

    # Panel D: Stress gene frequency
    ax = axes[1, 1]
    stress_only = amr_df[amr_df["Type"] == "STRESS"]
    stress_freq = stress_only.groupby("Element symbol")["source_file"].nunique().sort_values(ascending=True)
    ax.barh(range(len(stress_freq)), stress_freq.values,
            color=TYPE_COLORS["STRESS"], height=0.7)
    ax.set_yticks(range(len(stress_freq)))
    ax.set_yticklabels(stress_freq.index, fontsize=9)
    ax.set_xlabel("Number of Plasmids", fontsize=10)
    ax.set_title("D. Stress Gene Prevalence", fontsize=11, fontweight="bold")

    fig.suptitle(
        "AMR & Stress Gene Analysis \u2014 IncX Test Plasmids (n=22)",
        fontsize=14, fontweight="bold", y=1.01,
    )
    plt.tight_layout()
    return fig


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("pLIN Test \u2014 AMR Integration & Cladogram")
    print("=" * 70)

    # Part 1: Integration
    integrated, amr_df = integrate_plin_amr()

    # Part 2: Load sequences and compute clustering
    print(f"\n{'=' * 70}")
    print("PART 2: Generating AMR-annotated cladograms")
    print("=" * 70)

    print("\n  Loading sequences ...")
    records = load_sequences()
    print(f"  Loaded {len(records)} sequences")

    print("  Computing 4-mer vectors ...")
    vectors = compute_kmer_vectors(records, k=4)

    print("  Computing distances and clustering ...")
    dist_condensed = pdist(vectors, metric="cosine")
    Z = linkage(dist_condensed, method="single")

    plin_codes, cluster_assignments = get_plin_codes(records, Z)
    strain_clusters = list(cluster_assignments["F"])
    labels = [r["source_file"].replace(".fasta", "") for r in records]

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Figure 1: Cladogram with AMR heatmap
    print("  Generating cladogram with AMR profile ...")
    fig1 = plot_cladogram_amr(Z, labels, plin_codes, strain_clusters,
                               integrated, amr_df, records)
    for ext in ["png", "pdf"]:
        fig1.savefig(os.path.join(OUTPUT_DIR, f"test_cladogram_AMR.{ext}"),
                     dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig1)
    print("  Saved: test_cladogram_AMR.png/pdf")

    # Figure 2: AMR summary
    print("  Generating AMR summary figure ...")
    fig2 = plot_amr_summary(integrated, amr_df)
    for ext in ["png", "pdf"]:
        fig2.savefig(os.path.join(OUTPUT_DIR, f"test_AMR_summary.{ext}"),
                     dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig2)
    print("  Saved: test_AMR_summary.png/pdf")

    print(f"\n{'=' * 70}")
    print("All outputs saved to:", OUTPUT_DIR)
    print("  - test_pLIN_AMR_integrated.tsv")
    print("  - test_cladogram_AMR.png/pdf")
    print("  - test_AMR_summary.png/pdf")
    print("=" * 70)


if __name__ == "__main__":
    main()
