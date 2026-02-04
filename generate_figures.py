#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Generate publication-quality figures for the pLIN + AMRFinderPlus manuscript.
Produces 6 individual figures + 1 composite multi-panel figure.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap
from collections import Counter
import seaborn as sns

# ── Setup ─────────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
INTEGRATED = os.path.join(BASE_DIR, "output", "integrated", "pLIN_AMR_integrated.tsv")
PLIN_FILE = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")
AMR_RAW = os.path.join(BASE_DIR, "output", "amrfinder", "amrfinder_all_plasmids.tsv")
FIG_DIR = os.path.join(BASE_DIR, "output", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Style
sns.set_style("whitegrid")
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.2,
})

INC_COLORS = {"IncFII": "#2196F3", "IncN": "#FF9800", "IncX1": "#4CAF50"}
INC_ORDER = ["IncFII", "IncN", "IncX1"]

# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading data ...")
merged = pd.read_csv(INTEGRATED, sep="\t")
plin = pd.read_csv(PLIN_FILE, sep="\t")
amr_raw = pd.read_csv(AMR_RAW, sep="\t")

print(f"  Integrated: {len(merged)} plasmids")
print(f"  AMR raw: {len(amr_raw)} detections")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Dataset Overview — Inc type composition, size distribution, pLIN diversity
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating Figure 1: Dataset Overview ...")

fig1, axes1 = plt.subplots(1, 3, figsize=(14, 4.5))

# Panel A: Inc type composition (donut chart)
ax = axes1[0]
inc_counts = plin["inc_type"].value_counts().reindex(INC_ORDER)
wedges, texts, autotexts = ax.pie(
    inc_counts.values,
    labels=None,
    colors=[INC_COLORS[i] for i in INC_ORDER],
    autopct=lambda p: f"{p:.1f}%\n(n={int(p*sum(inc_counts)/100):,})",
    startangle=90,
    pctdistance=0.75,
    wedgeprops=dict(width=0.45, edgecolor="white", linewidth=2),
)
for at in autotexts:
    at.set_fontsize(8)
    at.set_fontweight("bold")
ax.legend(INC_ORDER, loc="lower center", ncol=3, fontsize=9, frameon=False,
          bbox_to_anchor=(0.5, -0.05))
ax.set_title("A. Dataset Composition", fontweight="bold", pad=15)

# Panel B: Plasmid size distribution by Inc type
ax = axes1[1]
for inc in INC_ORDER:
    sub = plin[plin["inc_type"] == inc]
    ax.hist(sub["length_bp"] / 1000, bins=50, alpha=0.6, color=INC_COLORS[inc],
            label=f"{inc} (n={len(sub):,})", edgecolor="white", linewidth=0.5)
ax.set_xlabel("Plasmid Length (kb)")
ax.set_ylabel("Count")
ax.set_title("B. Size Distribution", fontweight="bold")
ax.legend(frameon=True, fancybox=True, shadow=False, fontsize=8)
ax.set_xlim(0, 400)

# Panel C: pLIN diversity per Inc type (unique codes at each bin level)
ax = axes1[2]
bins = ["bin_A", "bin_B", "bin_C", "bin_D", "bin_E", "bin_F"]
bin_labels_short = ["A\nFamily", "B\nSubfamily", "C\nCluster", "D\nSubcluster", "E\nClone", "F\nStrain"]
x = np.arange(len(bins))
width = 0.25
for i, inc in enumerate(INC_ORDER):
    sub = plin[plin["inc_type"] == inc]
    n_clusters = [sub[b].nunique() for b in bins]
    bars = ax.bar(x + i * width, n_clusters, width, color=INC_COLORS[inc],
                  label=inc, edgecolor="white", linewidth=0.5)
    for bar, val in zip(bars, n_clusters):
        if val > 5:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 15,
                    str(val), ha="center", va="bottom", fontsize=6.5, fontweight="bold")
ax.set_xticks(x + width)
ax.set_xticklabels(bin_labels_short, fontsize=8)
ax.set_ylabel("Unique Clusters")
ax.set_title("C. Hierarchical pLIN Diversity", fontweight="bold")
ax.legend(frameon=True, fancybox=True, shadow=False, fontsize=8)
ax.set_yscale("log")
ax.set_ylim(0.8, 3000)

fig1.tight_layout(w_pad=3)
fig1.savefig(os.path.join(FIG_DIR, "Figure1_dataset_overview.png"))
fig1.savefig(os.path.join(FIG_DIR, "Figure1_dataset_overview.pdf"))
plt.close(fig1)
print("  Figure 1 saved.")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: AMR Gene Prevalence Overview
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 2: AMR Prevalence ...")

fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: Stacked bar — AMR/VIR/STRESS prevalence by Inc type
ax = axes2[0]
categories = ["AMR", "Virulence", "Stress"]
cat_colors = ["#E53935", "#8E24AA", "#FF8F00"]
bar_data = []
for inc in INC_ORDER:
    sub = merged[merged["inc_type_x"] == inc]
    n = len(sub)
    bar_data.append([
        (sub["n_amr_genes"] > 0).sum() / n * 100,
        (sub["n_vir_genes"] > 0).sum() / n * 100,
        (sub["n_stress_genes"] > 0).sum() / n * 100,
    ])
bar_data = np.array(bar_data)
x = np.arange(len(INC_ORDER))
width = 0.22
for j, (cat, col) in enumerate(zip(categories, cat_colors)):
    bars = ax.bar(x + j * width, bar_data[:, j], width, color=col, label=cat,
                  edgecolor="white", linewidth=0.5)
    for bar, val in zip(bars, bar_data[:, j]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f"{val:.0f}%", ha="center", va="bottom", fontsize=7, fontweight="bold")
ax.set_xticks(x + width)
ax.set_xticklabels(INC_ORDER, fontsize=10)
ax.set_ylabel("Prevalence (%)")
ax.set_ylim(0, 100)
ax.set_title("A. Gene Prevalence by Inc Type", fontweight="bold")
ax.legend(frameon=True, fancybox=True, shadow=False)

# Panel B: Top 15 AMR genes — horizontal bar chart
ax = axes2[1]
amr_pos = merged[merged["n_amr_genes"] > 0]
all_genes = []
for g in amr_pos["amr_genes"]:
    if g != "none":
        all_genes.extend(g.split("; "))
gene_counts = Counter(all_genes)
top15 = gene_counts.most_common(15)
genes_list = [g for g, _ in top15][::-1]
counts_list = [c for _, c in top15][::-1]
n_amr_pos = len(amr_pos)
pcts = [c / n_amr_pos * 100 for c in counts_list]

colors_bar = plt.cm.Reds(np.linspace(0.3, 0.85, len(genes_list)))
bars = ax.barh(range(len(genes_list)), pcts, color=colors_bar, edgecolor="white", linewidth=0.5)
ax.set_yticks(range(len(genes_list)))
ax.set_yticklabels([f"$\\it{{{g}}}$" for g in genes_list], fontsize=8)
ax.set_xlabel("% of AMR+ Plasmids")
ax.set_title("B. Top 15 AMR Genes", fontweight="bold")
for bar, pct, cnt in zip(bars, pcts, counts_list):
    ax.text(pct + 0.5, bar.get_y() + bar.get_height()/2,
            f"{pct:.1f}% (n={cnt:,})", va="center", fontsize=7)
ax.set_xlim(0, max(pcts) * 1.35)

# Panel C: Drug class distribution — horizontal bar chart
ax = axes2[2]
all_classes = []
for c in merged["amr_classes"]:
    if c != "none":
        all_classes.extend(c.split("; "))
class_counts = Counter(all_classes)
top12 = class_counts.most_common(12)
cls_list = [c for c, _ in top12][::-1]
cls_counts = [c for _, c in top12][::-1]

colors_cls = plt.cm.Blues(np.linspace(0.3, 0.85, len(cls_list)))
bars = ax.barh(range(len(cls_list)), cls_counts, color=colors_cls, edgecolor="white", linewidth=0.5)
ax.set_yticks(range(len(cls_list)))
ax.set_yticklabels(cls_list, fontsize=8)
ax.set_xlabel("Number of Plasmids")
ax.set_title("C. AMR Drug Classes", fontweight="bold")
for bar, cnt in zip(bars, cls_counts):
    ax.text(cnt + 20, bar.get_y() + bar.get_height()/2,
            f"n={cnt:,}", va="center", fontsize=7)
ax.set_xlim(0, max(cls_counts) * 1.2)

fig2.tight_layout(w_pad=3)
fig2.savefig(os.path.join(FIG_DIR, "Figure2_AMR_prevalence.png"))
fig2.savefig(os.path.join(FIG_DIR, "Figure2_AMR_prevalence.pdf"))
plt.close(fig2)
print("  Figure 2 saved.")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: Clinically Critical Resistance Genes
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 3: Critical Resistance Genes ...")

fig3, axes3 = plt.subplots(2, 2, figsize=(12, 10))

# Collect all AMR genes flat
all_genes_flat = []
for g in merged["amr_genes"]:
    if g != "none":
        all_genes_flat.extend(g.split("; "))

critical_data = {
    "Carbapenemases": {
        "patterns": ["blaKPC", "blaNDM", "blaOXA-48", "blaVIM", "blaIMP"],
        "color": "#D32F2F",
    },
    "ESBLs": {
        "patterns": ["blaCTX-M", "blaSHV", "blaTEM"],
        "color": "#F57C00",
    },
    "Colistin (mcr)": {
        "patterns": ["mcr-"],
        "color": "#7B1FA2",
    },
    "PMQR (Quinolone)": {
        "patterns": ["qnr", "aac(6')-Ib-cr", "oqxA", "oqxB"],
        "color": "#1976D2",
    },
}

for idx, (cat_name, cat_info) in enumerate(critical_data.items()):
    ax = axes3[idx // 2, idx % 2]
    matching = [g for g in all_genes_flat if any(p in g for p in cat_info["patterns"])]
    gene_counts_crit = Counter(matching)
    top10 = gene_counts_crit.most_common(10)
    if not top10:
        ax.text(0.5, 0.5, "None detected", ha="center", va="center", fontsize=14,
                transform=ax.transAxes)
        ax.set_title(cat_name, fontweight="bold")
        continue

    g_names = [g for g, _ in top10][::-1]
    g_counts = [c for _, c in top10][::-1]

    color_grad = plt.cm.Reds(np.linspace(0.25, 0.85, len(g_names))) if "Carbapenem" in cat_name else \
                 plt.cm.Oranges(np.linspace(0.25, 0.85, len(g_names))) if "ESBL" in cat_name else \
                 plt.cm.Purples(np.linspace(0.25, 0.85, len(g_names))) if "Colistin" in cat_name else \
                 plt.cm.Blues(np.linspace(0.25, 0.85, len(g_names)))

    bars = ax.barh(range(len(g_names)), g_counts, color=color_grad,
                   edgecolor="white", linewidth=0.5)
    ax.set_yticks(range(len(g_names)))
    ax.set_yticklabels([f"$\\it{{{g}}}$" for g in g_names], fontsize=9)
    ax.set_xlabel("Detections")
    ax.set_title(f"{cat_name} (n={sum(g_counts):,} total)", fontweight="bold", fontsize=11)
    for bar, cnt in zip(bars, g_counts):
        ax.text(cnt + max(g_counts)*0.02, bar.get_y() + bar.get_height()/2,
                f"{cnt:,}", va="center", fontsize=8)
    ax.set_xlim(0, max(g_counts) * 1.18)

fig3.suptitle("Clinically Critical Resistance Determinants", fontsize=14, fontweight="bold", y=1.01)
fig3.tight_layout()
fig3.savefig(os.path.join(FIG_DIR, "Figure3_critical_AMR.png"))
fig3.savefig(os.path.join(FIG_DIR, "Figure3_critical_AMR.pdf"))
plt.close(fig3)
print("  Figure 3 saved.")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: pLIN Lineages as AMR Vehicles — Heatmap
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 4: pLIN Lineage AMR Heatmap ...")

# Get top 10 pLIN lineages by AMR+ count
amr_pos_merged = merged[merged["n_amr_genes"] > 0]
plin_counts = amr_pos_merged["pLIN"].value_counts().head(10)
top_plins = plin_counts.index.tolist()

# Get top 15 AMR genes overall
top15_genes = [g for g, _ in Counter(all_genes).most_common(15)]

# Build heatmap matrix: prevalence of each gene in each lineage
heatmap_data = np.zeros((len(top_plins), len(top15_genes)))
lineage_info = []

for i, pcode in enumerate(top_plins):
    sub = merged[(merged["pLIN"] == pcode) & (merged["n_amr_genes"] > 0)]
    n = len(sub)
    inc_types = ", ".join(sorted(sub["inc_type_x"].unique()))
    lineage_info.append(f"{pcode}\n({inc_types}, n={n})")
    genes_in_lineage = []
    for g in sub["amr_genes"]:
        if g != "none":
            genes_in_lineage.extend(g.split("; "))
    gene_counter = Counter(genes_in_lineage)
    for j, gene in enumerate(top15_genes):
        heatmap_data[i, j] = gene_counter.get(gene, 0) / n * 100

fig4, ax4 = plt.subplots(figsize=(14, 7))
cmap = LinearSegmentedColormap.from_list("custom", ["#FFFFFF", "#FFCDD2", "#E53935", "#B71C1C"])
im = ax4.imshow(heatmap_data, cmap=cmap, aspect="auto", vmin=0, vmax=100)

ax4.set_xticks(range(len(top15_genes)))
ax4.set_xticklabels([f"$\\it{{{g}}}$" for g in top15_genes], rotation=45, ha="right", fontsize=9)
ax4.set_yticks(range(len(top_plins)))
ax4.set_yticklabels(lineage_info, fontsize=8)

# Annotate cells
for i in range(len(top_plins)):
    for j in range(len(top15_genes)):
        val = heatmap_data[i, j]
        if val > 0:
            color = "white" if val > 55 else "black"
            ax4.text(j, i, f"{val:.0f}", ha="center", va="center",
                     fontsize=7, color=color, fontweight="bold" if val > 50 else "normal")

cbar = plt.colorbar(im, ax=ax4, shrink=0.8, label="Prevalence (%)")
ax4.set_title("AMR Gene Prevalence (%) in Top pLIN Lineages", fontweight="bold", fontsize=13, pad=15)
ax4.set_xlabel("AMR Gene", fontsize=11)
ax4.set_ylabel("pLIN Lineage", fontsize=11)

fig4.tight_layout()
fig4.savefig(os.path.join(FIG_DIR, "Figure4_pLIN_AMR_heatmap.png"))
fig4.savefig(os.path.join(FIG_DIR, "Figure4_pLIN_AMR_heatmap.pdf"))
plt.close(fig4)
print("  Figure 4 saved.")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5: AMR Gene Burden Distribution + Virulence by Inc Type
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 5: AMR Burden & Virulence ...")

fig5, axes5 = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: Violin plot of AMR gene count per plasmid by Inc type
ax = axes5[0]
plot_data = []
for inc in INC_ORDER:
    sub = merged[(merged["inc_type_x"] == inc) & (merged["n_amr_genes"] > 0)]
    for val in sub["n_amr_genes"]:
        plot_data.append({"Inc Type": inc, "AMR Genes": val})
plot_df = pd.DataFrame(plot_data)

parts = ax.violinplot(
    [plot_df[plot_df["Inc Type"] == inc]["AMR Genes"].values for inc in INC_ORDER],
    positions=range(len(INC_ORDER)),
    showmeans=True, showmedians=True, showextrema=False,
)
for i, pc in enumerate(parts["bodies"]):
    pc.set_facecolor(INC_COLORS[INC_ORDER[i]])
    pc.set_alpha(0.7)
parts["cmeans"].set_color("black")
parts["cmedians"].set_color("red")

ax.set_xticks(range(len(INC_ORDER)))
ax.set_xticklabels(INC_ORDER)
ax.set_ylabel("AMR Genes per Plasmid")
ax.set_title("A. AMR Gene Burden\n(AMR+ plasmids only)", fontweight="bold")

# Add mean labels
for i, inc in enumerate(INC_ORDER):
    sub = merged[(merged["inc_type_x"] == inc) & (merged["n_amr_genes"] > 0)]
    mean_val = sub["n_amr_genes"].mean()
    ax.text(i, mean_val + 1, f"μ={mean_val:.1f}", ha="center", fontsize=8, fontweight="bold")

# Panel B: Virulence gene top 10 per Inc type (grouped bar)
ax = axes5[1]
vir_by_inc = {}
for inc in INC_ORDER:
    sub = merged[(merged["inc_type_x"] == inc) & (merged["n_vir_genes"] > 0)]
    vgenes = []
    for g in sub["vir_genes"]:
        if g != "none":
            vgenes.extend(g.split("; "))
    vir_by_inc[inc] = Counter(vgenes)

# Get top 8 virulence genes overall
all_vir = sum(vir_by_inc.values(), Counter())
top8_vir = [g for g, _ in all_vir.most_common(8)]

x = np.arange(len(top8_vir))
width = 0.25
for i, inc in enumerate(INC_ORDER):
    counts = [vir_by_inc[inc].get(g, 0) for g in top8_vir]
    ax.bar(x + i * width, counts, width, color=INC_COLORS[inc], label=inc,
           edgecolor="white", linewidth=0.5)
ax.set_xticks(x + width)
ax.set_xticklabels([f"$\\it{{{g}}}$" for g in top8_vir], rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Detections")
ax.set_title("B. Top Virulence Genes\nby Inc Type", fontweight="bold")
ax.legend(frameon=True, fancybox=True, shadow=False, fontsize=8)

# Panel C: Co-occurrence — AMR vs Virulence scatter
ax = axes5[2]
for inc in INC_ORDER:
    sub = merged[(merged["inc_type_x"] == inc) & (merged["n_total_hits"] > 0)]
    ax.scatter(sub["n_amr_genes"], sub["n_vir_genes"],
               c=INC_COLORS[inc], alpha=0.3, s=15, label=inc, edgecolors="none")
ax.set_xlabel("AMR Genes per Plasmid")
ax.set_ylabel("Virulence Genes per Plasmid")
ax.set_title("C. AMR–Virulence\nCo-occurrence", fontweight="bold")
ax.legend(frameon=True, fancybox=True, shadow=False, fontsize=8, markerscale=2)

fig5.tight_layout(w_pad=3)
fig5.savefig(os.path.join(FIG_DIR, "Figure5_AMR_burden_virulence.png"))
fig5.savefig(os.path.join(FIG_DIR, "Figure5_AMR_burden_virulence.pdf"))
plt.close(fig5)
print("  Figure 5 saved.")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 6: pLIN Hierarchical Structure Visualization
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 6: pLIN Hierarchical Structure ...")

fig6, axes6 = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Cluster count vs threshold (all Inc types combined + per type)
ax = axes6[0]
bins_list = ["bin_A", "bin_B", "bin_C", "bin_D", "bin_E", "bin_F"]
thresholds = [0.150, 0.100, 0.050, 0.020, 0.010, 0.001]
thresh_labels = ["A\n(0.150)", "B\n(0.100)", "C\n(0.050)", "D\n(0.020)", "E\n(0.010)", "F\n(0.001)"]

# Overall
overall_clusters = [plin[b].nunique() for b in bins_list]
ax.plot(range(len(bins_list)), overall_clusters, "ko-", linewidth=2.5, markersize=8,
        label=f"All (n={len(plin):,})", zorder=5)

for inc in INC_ORDER:
    sub = plin[plin["inc_type"] == inc]
    clusters = [sub[b].nunique() for b in bins_list]
    ax.plot(range(len(bins_list)), clusters, "o-", color=INC_COLORS[inc],
            linewidth=1.5, markersize=6, label=f"{inc} (n={len(sub):,})")

ax.set_xticks(range(len(bins_list)))
ax.set_xticklabels(thresh_labels, fontsize=9)
ax.set_ylabel("Number of Clusters")
ax.set_yscale("log")
ax.set_title("A. Clustering Resolution vs. Threshold", fontweight="bold")
ax.legend(frameon=True, fancybox=True, shadow=False, fontsize=8)
ax.grid(True, alpha=0.3)

# Panel B: Strain-level cluster size distribution
ax = axes6[1]
cluster_sizes = plin["pLIN"].value_counts()
size_dist = cluster_sizes.value_counts().sort_index()

ax.bar(size_dist.index[:30], size_dist.values[:30], color="#5C6BC0",
       edgecolor="white", linewidth=0.3)
ax.set_xlabel("Cluster Size (number of plasmids)")
ax.set_ylabel("Number of pLIN Codes")
ax.set_title("B. Strain-Level Cluster Size Distribution\n(Bin F, d ≤ 0.001)", fontweight="bold")

# Annotate singletons
n_sing = size_dist.get(1, 0)
ax.annotate(f"Singletons\nn={n_sing:,} ({n_sing/len(cluster_sizes)*100:.1f}%)",
            xy=(1, n_sing), xytext=(5, n_sing * 0.8),
            arrowprops=dict(arrowstyle="->", color="red"),
            fontsize=9, color="red", fontweight="bold")

fig6.tight_layout(w_pad=3)
fig6.savefig(os.path.join(FIG_DIR, "Figure6_pLIN_hierarchy.png"))
fig6.savefig(os.path.join(FIG_DIR, "Figure6_pLIN_hierarchy.pdf"))
plt.close(fig6)
print("  Figure 6 saved.")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 7: Composite Multi-Panel Figure for Main Manuscript
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 7: Composite Manuscript Figure ...")

fig7 = plt.figure(figsize=(18, 20))
gs = gridspec.GridSpec(4, 3, figure=fig7, hspace=0.35, wspace=0.35,
                       height_ratios=[1, 1, 1.2, 1])

# ── Row 1: Dataset overview ──
# A: Inc composition donut
ax = fig7.add_subplot(gs[0, 0])
wedges, texts, autotexts = ax.pie(
    inc_counts.values, labels=None,
    colors=[INC_COLORS[i] for i in INC_ORDER],
    autopct=lambda p: f"{p:.1f}%",
    startangle=90, pctdistance=0.75,
    wedgeprops=dict(width=0.45, edgecolor="white", linewidth=2),
)
for at in autotexts:
    at.set_fontsize(8)
ax.legend(INC_ORDER, loc="lower center", ncol=3, fontsize=7, frameon=False,
          bbox_to_anchor=(0.5, -0.08))
ax.set_title("A. Dataset Composition\n(n=6,346)", fontweight="bold", fontsize=10)

# B: Size distribution
ax = fig7.add_subplot(gs[0, 1])
for inc in INC_ORDER:
    sub = plin[plin["inc_type"] == inc]
    ax.hist(sub["length_bp"] / 1000, bins=50, alpha=0.6, color=INC_COLORS[inc],
            label=inc, edgecolor="white", linewidth=0.3)
ax.set_xlabel("Length (kb)", fontsize=9)
ax.set_ylabel("Count", fontsize=9)
ax.set_title("B. Plasmid Size Distribution", fontweight="bold", fontsize=10)
ax.legend(frameon=True, fontsize=7)
ax.set_xlim(0, 400)

# C: Clustering resolution
ax = fig7.add_subplot(gs[0, 2])
ax.plot(range(len(bins_list)), overall_clusters, "ko-", linewidth=2, markersize=6, label="All", zorder=5)
for inc in INC_ORDER:
    sub = plin[plin["inc_type"] == inc]
    clusters = [sub[b].nunique() for b in bins_list]
    ax.plot(range(len(bins_list)), clusters, "o-", color=INC_COLORS[inc],
            linewidth=1.2, markersize=4, label=inc)
ax.set_xticks(range(len(bins_list)))
ax.set_xticklabels(["A", "B", "C", "D", "E", "F"], fontsize=8)
ax.set_ylabel("Clusters", fontsize=9)
ax.set_yscale("log")
ax.set_title("C. Hierarchical Resolution", fontweight="bold", fontsize=10)
ax.legend(frameon=True, fontsize=7)

# ── Row 2: AMR prevalence ──
# D: AMR/VIR/STRESS prevalence
ax = fig7.add_subplot(gs[1, 0])
categories_short = ["AMR", "VIR", "Stress"]
x = np.arange(len(INC_ORDER))
width = 0.22
for j, (cat, col) in enumerate(zip(categories_short, cat_colors)):
    vals = bar_data[:, j]
    ax.bar(x + j * width, vals, width, color=col, label=cat, edgecolor="white", linewidth=0.3)
    for xi, v in zip(x + j * width, vals):
        ax.text(xi, v + 1.5, f"{v:.0f}%", ha="center", fontsize=6, fontweight="bold")
ax.set_xticks(x + width)
ax.set_xticklabels(INC_ORDER, fontsize=9)
ax.set_ylabel("Prevalence (%)", fontsize=9)
ax.set_ylim(0, 100)
ax.set_title("D. Gene Prevalence by Inc Type", fontweight="bold", fontsize=10)
ax.legend(frameon=True, fontsize=7)

# E: Top 10 AMR genes
ax = fig7.add_subplot(gs[1, 1])
top10_genes = gene_counts.most_common(10)
g10_names = [g for g, _ in top10_genes][::-1]
g10_counts = [c for _, c in top10_genes][::-1]
g10_pcts = [c / n_amr_pos * 100 for c in g10_counts]
colors10 = plt.cm.Reds(np.linspace(0.3, 0.85, len(g10_names)))
bars = ax.barh(range(len(g10_names)), g10_pcts, color=colors10, edgecolor="white", linewidth=0.3)
ax.set_yticks(range(len(g10_names)))
ax.set_yticklabels([f"$\\it{{{g}}}$" for g in g10_names], fontsize=7)
ax.set_xlabel("% of AMR+ plasmids", fontsize=9)
ax.set_title("E. Top 10 AMR Genes", fontweight="bold", fontsize=10)
for bar, pct in zip(bars, g10_pcts):
    ax.text(pct + 0.3, bar.get_y() + bar.get_height()/2, f"{pct:.1f}%", va="center", fontsize=6)

# F: Drug classes
ax = fig7.add_subplot(gs[1, 2])
top8_cls = class_counts.most_common(8)
c8_names = [c for c, _ in top8_cls][::-1]
c8_counts = [c for _, c in top8_cls][::-1]
colors8 = plt.cm.Blues(np.linspace(0.3, 0.85, len(c8_names)))
bars = ax.barh(range(len(c8_names)), c8_counts, color=colors8, edgecolor="white", linewidth=0.3)
ax.set_yticks(range(len(c8_names)))
ax.set_yticklabels(c8_names, fontsize=7)
ax.set_xlabel("Plasmids", fontsize=9)
ax.set_title("F. AMR Drug Classes", fontweight="bold", fontsize=10)

# ── Row 3: pLIN-AMR heatmap (spans full width) ──
ax = fig7.add_subplot(gs[2, :])
im = ax.imshow(heatmap_data, cmap=cmap, aspect="auto", vmin=0, vmax=100)
ax.set_xticks(range(len(top15_genes)))
ax.set_xticklabels([f"$\\it{{{g}}}$" for g in top15_genes], rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(len(top_plins)))
short_labels = []
for pcode in top_plins:
    sub = merged[(merged["pLIN"] == pcode) & (merged["n_amr_genes"] > 0)]
    incs = ",".join(sorted(sub["inc_type_x"].unique()))
    short_labels.append(f"{pcode} ({incs}, n={len(sub)})")
ax.set_yticklabels(short_labels, fontsize=7)
for i in range(len(top_plins)):
    for j in range(len(top15_genes)):
        val = heatmap_data[i, j]
        if val > 0:
            color = "white" if val > 55 else "black"
            ax.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=6.5,
                    color=color, fontweight="bold" if val > 50 else "normal")
cbar = plt.colorbar(im, ax=ax, shrink=0.6, label="Prevalence (%)", pad=0.02)
cbar.ax.tick_params(labelsize=7)
ax.set_title("G. AMR Gene Prevalence (%) in Top pLIN Lineages", fontweight="bold", fontsize=11)

# ── Row 4: Critical genes + co-occurrence ──
# H: Carbapenemases
ax = fig7.add_subplot(gs[3, 0])
carb_genes = [g for g in all_genes_flat if any(p in g for p in ["blaKPC", "blaNDM", "blaOXA-48", "blaVIM", "blaIMP"])]
carb_counts = Counter(carb_genes)
top8_carb = carb_counts.most_common(8)
if top8_carb:
    cn = [g for g, _ in top8_carb][::-1]
    cv = [c for _, c in top8_carb][::-1]
    colors_c = plt.cm.Reds(np.linspace(0.3, 0.85, len(cn)))
    bars = ax.barh(range(len(cn)), cv, color=colors_c, edgecolor="white", linewidth=0.3)
    ax.set_yticks(range(len(cn)))
    ax.set_yticklabels([f"$\\it{{{g}}}$" for g in cn], fontsize=7)
    ax.set_xlabel("Detections", fontsize=9)
    for bar, cnt in zip(bars, cv):
        ax.text(cnt + 5, bar.get_y() + bar.get_height()/2, str(cnt), va="center", fontsize=7)
ax.set_title("H. Carbapenemases", fontweight="bold", fontsize=10)

# I: Colistin + PMQR
ax = fig7.add_subplot(gs[3, 1])
mcr_genes = [g for g in all_genes_flat if "mcr-" in g]
mcr_counts = Counter(mcr_genes)
top6_mcr = mcr_counts.most_common(6)
if top6_mcr:
    mn = [g for g, _ in top6_mcr][::-1]
    mv = [c for _, c in top6_mcr][::-1]
    colors_m = plt.cm.Purples(np.linspace(0.3, 0.85, len(mn)))
    bars = ax.barh(range(len(mn)), mv, color=colors_m, edgecolor="white", linewidth=0.3)
    ax.set_yticks(range(len(mn)))
    ax.set_yticklabels([f"$\\it{{{g}}}$" for g in mn], fontsize=7)
    ax.set_xlabel("Detections", fontsize=9)
    for bar, cnt in zip(bars, mv):
        ax.text(cnt + 1, bar.get_y() + bar.get_height()/2, str(cnt), va="center", fontsize=7)
ax.set_title("I. Colistin Resistance (mcr)", fontweight="bold", fontsize=10)

# J: AMR-Virulence co-occurrence scatter
ax = fig7.add_subplot(gs[3, 2])
for inc in INC_ORDER:
    sub = merged[(merged["inc_type_x"] == inc) & (merged["n_total_hits"] > 0)]
    ax.scatter(sub["n_amr_genes"], sub["n_vir_genes"],
               c=INC_COLORS[inc], alpha=0.25, s=10, label=inc, edgecolors="none")
ax.set_xlabel("AMR Genes", fontsize=9)
ax.set_ylabel("Virulence Genes", fontsize=9)
ax.set_title("J. AMR–Virulence Co-occurrence", fontweight="bold", fontsize=10)
ax.legend(frameon=True, fontsize=7, markerscale=2)

fig7.savefig(os.path.join(FIG_DIR, "Figure7_composite_manuscript.png"))
fig7.savefig(os.path.join(FIG_DIR, "Figure7_composite_manuscript.pdf"))
plt.close(fig7)
print("  Figure 7 saved.")


# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print(f"All figures saved to: {FIG_DIR}")
print(f"{'='*60}")
for f in sorted(os.listdir(FIG_DIR)):
    fpath = os.path.join(FIG_DIR, f)
    size_kb = os.path.getsize(fpath) / 1024
    print(f"  {f:<45s} {size_kb:>8.1f} KB")
print(f"{'='*60}")
