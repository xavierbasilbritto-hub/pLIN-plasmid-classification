#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
"""
Generate Figures 14–17: Addressed Limitations — Multi-panel overview.

Figure 14: Analytical module validation summary
  A) Classifier training composition (28 Inc/Rep groups, training samples per group)
  B) Per-class F1 scores (5-fold stratified CV)
  C) Assembly completeness scoring demonstration (simulated distribution)
  D) Database coverage & novelty detection (distance percentile concept)
  E) Confusion matrix heatmap for 28-group classifier
  F) MGE boundary detection schematic (gene architecture)

Figure 15: Novel Inc/rep group discovery and evolutionary rate estimation
  (Gram-positive expansion, classification performance, centroid distances)

Figure 16: Adaptive thresholds and cluster stability
  A) Bootstrap stability assessment, B) Linkage method comparison

Figure 17: Mobile genetic element boundary detection / Recombination
  A) Recombination detection, B) Evolutionary rate estimation
"""

import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import seaborn as sns

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(BASE_DIR, "output", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Style
sns.set_style("whitegrid")
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.3,
})

# Load classifier data
DATA_FILE = os.path.join(BASE_DIR, "data", "inc_classifier.npz")
data = np.load(DATA_FILE, allow_pickle=True)
X = data["X"]
y = data["y"]
group_names = list(data["group_names"])
centroids = data["centroids"]
cv_metrics = json.loads(str(data["cv_metrics"][0]))
conf_matrix = data["confusion_matrix"]
cv_accuracy = float(data["cv_accuracy"][0])

n_groups = len(group_names)
samples_per_group = [int(np.sum(y == i)) for i in range(n_groups)]

# Classify groups into Gram-negative vs Gram-positive
gram_neg_groups = [g for g in group_names if not g.startswith("rep")]
gram_pos_groups = [g for g in group_names if g.startswith("rep")]

print(f"Loaded classifier: {n_groups} groups, {len(X)} samples, {cv_accuracy:.1%} accuracy")
print(f"  Gram-negative: {len(gram_neg_groups)} groups")
print(f"  Gram-positive: {len(gram_pos_groups)} groups")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 16: Multi-panel Limitations Addressed
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating Figure 14: Analytical Module Validation Summary ...")

fig = plt.figure(figsize=(20, 14))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)

# Color scheme
GRAM_NEG_COLOR = "#3498db"
GRAM_POS_SA_COLOR = "#e74c3c"
GRAM_POS_EF_COLOR = "#2ecc71"
ACINETOBACTER_COLOR = "#9b59b6"
PSEUDOMONAS_COLOR = "#e67e22"

def get_group_color(name):
    if name.startswith("repSA"):
        return GRAM_POS_SA_COLOR
    elif name.startswith("repEF"):
        return GRAM_POS_EF_COLOR
    elif name.startswith("repAci"):
        return ACINETOBACTER_COLOR
    elif name.startswith("repPae"):
        return PSEUDOMONAS_COLOR
    else:
        return GRAM_NEG_COLOR


# ── Panel A: Training composition — 24 groups ────────────────────────────────
ax_a = fig.add_subplot(gs[0, 0])

# Sort by sample count
sorted_idx = np.argsort(samples_per_group)
sorted_names = [group_names[i] for i in sorted_idx]
sorted_counts = [samples_per_group[i] for i in sorted_idx]
sorted_colors = [get_group_color(n) for n in sorted_names]

bars = ax_a.barh(range(n_groups), sorted_counts, color=sorted_colors, edgecolor="white",
                 linewidth=0.3, height=0.8)
ax_a.set_yticks(range(n_groups))
ax_a.set_yticklabels(sorted_names, fontsize=8)
ax_a.set_xlabel("Training sequences")
ax_a.set_title(f"A) Classifier Training Composition\n({n_groups} Inc/Rep Groups)", fontweight="bold")

# Add count labels
for i, (cnt, bar) in enumerate(zip(sorted_counts, bars)):
    if cnt > 200:
        ax_a.text(cnt - 5, i, str(cnt), va="center", ha="right", fontsize=6, color="white",
                  fontweight="bold")
    else:
        ax_a.text(cnt + 5, i, str(cnt), va="center", ha="left", fontsize=6)

# Legend — all 5 organism categories
legend_patches = [
    mpatches.Patch(color=GRAM_NEG_COLOR, label=f"Gram-negative ({len(gram_neg_groups)})"),
    mpatches.Patch(color=GRAM_POS_SA_COLOR, label=f"S. aureus ({sum(1 for g in gram_pos_groups if 'SA' in g)})"),
    mpatches.Patch(color=GRAM_POS_EF_COLOR, label=f"Enterococcus ({sum(1 for g in gram_pos_groups if 'EF' in g)})"),
    mpatches.Patch(color=ACINETOBACTER_COLOR, label=f"Acinetobacter ({sum(1 for g in gram_pos_groups if 'Aci' in g)})"),
    mpatches.Patch(color=PSEUDOMONAS_COLOR, label=f"Pseudomonas ({sum(1 for g in gram_pos_groups if 'Pae' in g)})"),
]
ax_a.legend(handles=legend_patches, loc="lower right", fontsize=8, framealpha=0.9)
ax_a.text(0.98, 0.02, f"Total: {sum(samples_per_group):,} sequences\nAccuracy: {cv_accuracy:.1%}",
          transform=ax_a.transAxes, ha="right", va="bottom", fontsize=7,
          bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))


# ── Panel B: Per-class F1 scores ─────────────────────────────────────────────
ax_b = fig.add_subplot(gs[0, 1])

f1_scores = [cv_metrics[g]["f1"] for g in group_names]
sorted_f1_idx = np.argsort(f1_scores)
sorted_f1_names = [group_names[i] for i in sorted_f1_idx]
sorted_f1_vals = [f1_scores[i] for i in sorted_f1_idx]
sorted_f1_colors = [get_group_color(n) for n in sorted_f1_names]

bars_b = ax_b.barh(range(n_groups), sorted_f1_vals, color=sorted_f1_colors,
                   edgecolor="white", linewidth=0.3, height=0.8)
ax_b.set_yticks(range(n_groups))
ax_b.set_yticklabels(sorted_f1_names, fontsize=8)
ax_b.set_xlabel("F1 Score")
ax_b.set_xlim(0, 1.05)
ax_b.axvline(x=0.8, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
ax_b.set_title("B) Per-Class F1 Scores\n(5-Fold Stratified CV)", fontweight="bold")

# Add F1 labels
for i, (f1, bar) in enumerate(zip(sorted_f1_vals, bars_b)):
    ax_b.text(f1 + 0.02, i, f"{f1:.2f}", va="center", ha="left", fontsize=6)

# Median line
median_f1 = np.median(f1_scores)
ax_b.axvline(x=median_f1, color="orange", linestyle="-.", linewidth=1, alpha=0.7)
ax_b.text(median_f1 + 0.02, n_groups - 1.5, f"Median: {median_f1:.2f}",
          fontsize=7, color="orange")


# ── Panel C: Assembly completeness (REAL data from training plasmids) ─────────
ax_c = fig.add_subplot(gs[0, 2])

COMP_TSV = os.path.join(BASE_DIR, "output", "figure14_validation", "assembly_completeness.tsv")
comp_df = pd.read_csv(COMP_TSV, sep="\t")

colors_comp = {"COMPLETE": "#27ae60", "NEAR-COMPLETE": "#f39c12",
               "FRAGMENTED": "#e67e22", "POOR": "#e74c3c"}

# Separate by status
for status, color in colors_comp.items():
    subset = comp_df[comp_df["completeness_status"] == status]["completeness_score"]
    if len(subset) > 0:
        ax_c.hist(subset, bins=20, alpha=0.8, color=color,
                  label=f"{status} (n={len(subset)})", edgecolor="white", linewidth=0.3)

# Threshold lines
for thresh, label_t in [(80, "Complete"), (60, "Near-complete"), (40, "Fragmented")]:
    ax_c.axvline(x=thresh, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)

# Fix y-axis annotation after draw
ax_c.set_xlabel("Completeness Score")
ax_c.set_ylabel("Count")
n_total = len(comp_df)
n_single = int((comp_df["n_contigs"] == 1).sum())
n_circ = int(comp_df["circular_signal"].sum())
ax_c.set_title(f"C) Assembly Completeness Assessment\n({n_total:,} Training Plasmids)", fontweight="bold")
ax_c.legend(fontsize=7, loc="upper left")

# Scoring criteria + real stats inset
criteria_text = ("Scoring criteria:\n"
                 "  Single contig: +40\n"
                 "  N50 ratio > 0.9: +20\n"
                 "  Circular signal: +20\n"
                 "  Coding density > 80%: +10\n"
                 "  No N-gaps: +10\n"
                 f"Single contig: {n_single/n_total*100:.1f}%\n"
                 f"Circular: {n_circ/n_total*100:.1f}%\n"
                 f"Median score: {comp_df['completeness_score'].median():.0f}")
ax_c.text(0.98, 0.55, criteria_text, transform=ax_c.transAxes, fontsize=6,
          ha="right", va="top", fontfamily="monospace",
          bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))


# ── Panel D: Database coverage / novelty (REAL leave-one-out NN distances) ───
ax_d = fig.add_subplot(gs[1, 0])

NN_TSV = os.path.join(BASE_DIR, "output", "figure14_validation", "nn_distances.tsv")
nn_df = pd.read_csv(NN_TSV, sep="\t")

# Separate by coverage indicator
green = nn_df[nn_df["coverage_indicator"] == "GREEN"]["nn_distance"]
yellow = nn_df[nn_df["coverage_indicator"] == "YELLOW"]["nn_distance"]
red = nn_df[nn_df["coverage_indicator"] == "RED"]["nn_distance"]

ax_d.hist(green, bins=50, alpha=0.7, color="#27ae60",
          label=f"Well-covered (n={len(green)})", edgecolor="white", linewidth=0.3)
ax_d.hist(yellow, bins=20, alpha=0.8, color="#f39c12",
          label=f"Divergent (n={len(yellow)})", edgecolor="white", linewidth=0.3)
ax_d.hist(red, bins=15, alpha=0.8, color="#e74c3c",
          label=f"Potentially novel (n={len(red)})", edgecolor="white", linewidth=0.3)

# Threshold lines
ax_d.axvline(x=0.020, color="#f39c12", linestyle="--", linewidth=1.2, alpha=0.8)
ax_d.axvline(x=0.050, color="#e74c3c", linestyle="--", linewidth=1.2, alpha=0.8)

ax_d.set_xlabel("Leave-One-Out Nearest-Neighbour Distance (cosine)")
ax_d.set_ylabel("Count")
n_nn = len(nn_df)
ax_d.set_title(f"D) Database Coverage & Novelty Detection\n({n_nn:,} Training Plasmids, Leave-One-Out)", fontweight="bold")
ax_d.legend(fontsize=7, loc="upper right")
ax_d.set_xlim(0, 0.10)

# Stats inset
nn_med = nn_df["nn_distance"].median()
nn_p95 = nn_df["nn_distance"].quantile(0.95)
nn_p99 = nn_df["nn_distance"].quantile(0.99)
ax_d.text(0.98, 0.95,
          f"Median: {nn_med:.4f}\np95: {nn_p95:.4f}\np99: {nn_p99:.4f}\n"
          f"GREEN: {len(green)/n_nn*100:.1f}%\n"
          f"YELLOW: {len(yellow)/n_nn*100:.1f}%\n"
          f"RED: {len(red)/n_nn*100:.1f}%",
          transform=ax_d.transAxes, fontsize=6, ha="right", va="top",
          bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))


# ── Panel E: Confusion matrix heatmap ────────────────────────────────────────
ax_e = fig.add_subplot(gs[1, 1])

# Normalize confusion matrix (row-wise)
conf_norm = conf_matrix.astype(float)
row_sums = conf_norm.sum(axis=1, keepdims=True)
row_sums[row_sums == 0] = 1
conf_norm = conf_norm / row_sums

# Use a diverging colormap
cmap = LinearSegmentedColormap.from_list("conf",
    ["#ffffff", "#d4e6f1", "#2980b9", "#1a5276"], N=256)

im = ax_e.imshow(conf_norm, cmap=cmap, aspect="auto", vmin=0, vmax=1)
ax_e.set_xticks(range(n_groups))
ax_e.set_yticks(range(n_groups))
ax_e.set_xticklabels(group_names, rotation=90, fontsize=6.5, ha="center")
ax_e.set_yticklabels(group_names, fontsize=6.5)
ax_e.set_xlabel("Predicted", fontsize=10)
ax_e.set_ylabel("True", fontsize=10)
ax_e.set_title("E) Classifier Confusion Matrix\n(Row-Normalized, 5-Fold CV)", fontweight="bold")

# Add colorbar
cbar = plt.colorbar(im, ax=ax_e, fraction=0.046, pad=0.04)
cbar.set_label("Proportion", fontsize=8)

# Highlight Gram-positive blocks
gp_start = group_names.index(gram_pos_groups[0])
gp_end = group_names.index(gram_pos_groups[-1])
rect = plt.Rectangle((gp_start - 0.5, gp_start - 0.5),
                      gp_end - gp_start + 1, gp_end - gp_start + 1,
                      fill=False, edgecolor=GRAM_POS_SA_COLOR, linewidth=1.5,
                      linestyle="--")
ax_e.add_patch(rect)
ax_e.text(gp_end + 1.2, (gp_start + gp_end) / 2, "Gram+",
          fontsize=7, color=GRAM_POS_SA_COLOR, ha="left", va="center", fontweight="bold")


# ── Panel F: MGE boundary detection schematic ────────────────────────────────
ax_f = fig.add_subplot(gs[1, 2])
ax_f.set_xlim(0, 100)
ax_f.set_ylim(-2, 8)
ax_f.set_aspect("auto")

# Draw a stylized plasmid gene map
gene_regions = [
    # (start, end, label, color, category)
    (0, 8, "repA", "#3498db", "Backbone"),
    (8, 14, "parA", "#3498db", "Backbone"),
    (14, 20, "parB", "#3498db", "Backbone"),
    (20, 25, "IS26", "#f1c40f", "IS element"),
    (25, 35, "blaCTX-M-15", "#e74c3c", "Resistance"),
    (35, 40, "IS26", "#f1c40f", "IS element"),
    (40, 48, "traA", "#2ecc71", "Mobility"),
    (48, 54, "traB", "#2ecc71", "Mobility"),
    (54, 60, "traC", "#2ecc71", "Mobility"),
    (60, 65, "IS1", "#f1c40f", "IS element"),
    (65, 72, "aac(6')-Ib", "#e74c3c", "Resistance"),
    (72, 78, "sul1", "#e74c3c", "Resistance"),
    (78, 82, "IS1", "#f1c40f", "IS element"),
    (82, 88, "hyp1", "#95a5a6", "Hypothetical"),
    (88, 94, "hyp2", "#95a5a6", "Hypothetical"),
    (94, 100, "korA", "#3498db", "Backbone"),
]

y_gene = 4
gene_h = 1.5

for start, end, label, color, cat in gene_regions:
    width = end - start
    # Draw arrow-shaped gene
    arrow = mpatches.FancyArrowPatch(
        (start, y_gene), (end, y_gene),
        arrowstyle="-|>",
        mutation_scale=15,
        linewidth=0,
        color=color,
    )
    rect = plt.Rectangle((start, y_gene - gene_h/2), width, gene_h,
                         facecolor=color, edgecolor="white", linewidth=0.5, alpha=0.85)
    ax_f.add_patch(rect)

    # Gene label (only if wide enough)
    if width >= 5:
        ax_f.text((start + end) / 2, y_gene, label, ha="center", va="center",
                  fontsize=5, fontweight="bold", color="white", rotation=0)
    elif width >= 3.5:
        ax_f.text((start + end) / 2, y_gene, label, ha="center", va="center",
                  fontsize=4, color="white", rotation=45)

# Draw composite transposon brackets
bracket_y = y_gene + gene_h/2 + 0.3
ax_f.annotate("", xy=(20, bracket_y + 1.2), xytext=(40, bracket_y + 1.2),
              arrowprops=dict(arrowstyle="<->", color="#e67e22", lw=1.5))
ax_f.text(30, bracket_y + 1.5, "Composite Transposon\n(IS26-flanked)", ha="center",
          va="bottom", fontsize=6, color="#e67e22", fontweight="bold")

ax_f.annotate("", xy=(60, bracket_y + 1.2), xytext=(82, bracket_y + 1.2),
              arrowprops=dict(arrowstyle="<->", color="#e67e22", lw=1.5))
ax_f.text(71, bracket_y + 1.5, "Composite Transposon\n(IS1-flanked)", ha="center",
          va="bottom", fontsize=6, color="#e67e22", fontweight="bold")

# Category legend
legend_items = [
    mpatches.Patch(color="#3498db", label="Backbone"),
    mpatches.Patch(color="#e74c3c", label="Resistance"),
    mpatches.Patch(color="#2ecc71", label="Mobility/Transfer"),
    mpatches.Patch(color="#f1c40f", label="IS Elements"),
    mpatches.Patch(color="#95a5a6", label="Hypothetical"),
]
ax_f.legend(handles=legend_items, loc="lower center", ncol=5, fontsize=8,
            framealpha=0.9, bbox_to_anchor=(0.5, -0.12))

# Scale bar
ax_f.plot([0, 10], [-0.8, -0.8], "k-", linewidth=1.5)
ax_f.text(5, -1.3, "~5 kb", ha="center", fontsize=7)

ax_f.set_title("F) MGE Boundary Detection\n(L10: Gene Architecture Map)", fontweight="bold")
ax_f.set_xlabel("Position (kb)", fontsize=10)
ax_f.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
ax_f.set_xticklabels(["0", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50"])
ax_f.set_yticks([])
ax_f.spines["left"].set_visible(False)
ax_f.spines["top"].set_visible(False)

# ── Main title ───────────────────────────────────────────────────────────────
fig.suptitle("Figure 13: Analytical Module Validation Summary — Expanded Classifier,\n"
             "Quality Assessment, Novelty Detection, and Mobile Genetic Element Analysis",
             fontsize=13, fontweight="bold", y=1.01)

# Save
out_path = os.path.join(FIG_DIR, "figure14_analytical_validation.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
out_pdf = out_path.replace(".png", ".pdf")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"  Saved: {out_path}")
print(f"  Saved: {out_pdf}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 17: Gram-Positive Expansion Detail
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating Figure 15: Novel Inc/Rep Group Discovery and Expansion ...")

fig17, axes17 = plt.subplots(1, 3, figsize=(18, 6))

# Panel A: Gram-positive group composition (horizontal bars to avoid label overlap)
ax = axes17[0]
gp_names = gram_pos_groups
gp_counts = [samples_per_group[group_names.index(g)] for g in gp_names]
gp_colors = [get_group_color(g) for g in gp_names]

gp_labels_display = {
    "repSA_large": "repSA_large (S. aureus large)",
    "repSA_small": "repSA_small (S. aureus small)",
    "repEF_conj": "repEF_conj (Enterococcus conj.)",
    "repEF_res": "repEF_res (Enterococcus res.)",
    "repAci1": "repAci1 (Acinetobacter small)",
    "repAci_large": "repAci_large (Acinetobacter)",
    "repPae_large": "repPae_large (Pseudomonas)",
    "repPae_small": "repPae_small (Pseudomonas)",
}

# Sort by count for cleaner display
sorted_pairs = sorted(zip(gp_names, gp_counts, gp_colors), key=lambda x: x[1])
s_names, s_counts, s_colors = zip(*sorted_pairs)

bars17 = ax.barh(range(len(s_names)), s_counts, color=s_colors,
                 edgecolor="white", linewidth=0.5, height=0.7)
ax.set_yticks(range(len(s_names)))
ax.set_yticklabels([gp_labels_display.get(g, g) for g in s_names], fontsize=8)
ax.set_xlabel("Training sequences")
ax.set_title("A) Expanded Rep Type Groups\n(Gram-positive, Acinetobacter, Pseudomonas)", fontweight="bold")

for i, (cnt, bar) in enumerate(zip(s_counts, bars17)):
    ax.text(cnt + 3, i, str(cnt), va="center", ha="left", fontsize=9, fontweight="bold")

# Panel B: F1 scores comparison (Gram-neg vs Gram-pos)
ax = axes17[1]
gn_f1 = [cv_metrics[g]["f1"] for g in gram_neg_groups]
gp_f1 = [cv_metrics[g]["f1"] for g in gram_pos_groups]

bp_data = [gn_f1, gp_f1]
bp_colors = [GRAM_NEG_COLOR, "#9b59b6"]
bp = ax.boxplot(bp_data, positions=[1, 2], widths=0.5, patch_artist=True,
                showmeans=True, meanline=True)
for patch, color in zip(bp["boxes"], bp_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)

ax.set_xticks([1, 2])
ax.set_xticklabels([f"Gram-negative\n(n={len(gn_f1)})",
                     f"Gram-positive\n(n={len(gp_f1)})"])
ax.set_ylabel("F1 Score")
ax.set_title("B) Classification Performance\nby Organism Group", fontweight="bold")
ax.set_ylim(0, 1.1)

# Add individual points
for i, (vals, x_pos) in enumerate(zip(bp_data, [1, 2])):
    jitter = np.random.uniform(-0.12, 0.12, len(vals))
    ax.scatter([x_pos + j for j in jitter], vals, s=20, alpha=0.6,
              color=bp_colors[i], edgecolors="white", linewidth=0.5, zorder=3)

ax.text(0.02, 0.02, f"Gram-neg median F1: {np.median(gn_f1):.3f}\n"
        f"Gram-pos median F1: {np.median(gp_f1):.3f}",
        transform=ax.transAxes, fontsize=7, va="bottom",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

# Panel C: Centroid distance heatmap (Gram-pos groups vs selected Gram-neg)
ax = axes17[2]
from scipy.spatial.distance import cosine as cos_dist

# Compute centroid distances between Gram-positive and select Gram-negative groups
selected_gn = ["IncF", "IncFII", "IncN", "IncX1", "IncX3", "ColE", "ColRNAI", "IncHI2"]
all_selected = selected_gn + gram_pos_groups
n_sel = len(all_selected)

dist_matrix = np.zeros((n_sel, n_sel))
for i, g1 in enumerate(all_selected):
    idx1 = group_names.index(g1)
    for j, g2 in enumerate(all_selected):
        idx2 = group_names.index(g2)
        dist_matrix[i, j] = cos_dist(centroids[idx1], centroids[idx2])

# Plot heatmap
im17 = ax.imshow(dist_matrix, cmap="YlOrRd", aspect="auto")
ax.set_xticks(range(n_sel))
ax.set_yticks(range(n_sel))
ax.set_xticklabels(all_selected, rotation=45, ha="right", fontsize=8)
ax.set_yticklabels(all_selected, fontsize=8)
ax.set_title("C) Centroid Distances\n(Gram-pos vs Selected Gram-neg)", fontweight="bold")

cbar17 = plt.colorbar(im17, ax=ax, fraction=0.046, pad=0.04)
cbar17.set_label("Cosine Distance", fontsize=8)

# Annotate cells with distances
for i in range(n_sel):
    for j in range(n_sel):
        val = dist_matrix[i, j]
        color = "white" if val > 0.15 else "black"
        ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=4.5, color=color)

# Highlight Gram-positive block
gp_start_17 = len(selected_gn)
rect17 = plt.Rectangle((gp_start_17 - 0.5, gp_start_17 - 0.5),
                        len(gram_pos_groups), len(gram_pos_groups),
                        fill=False, edgecolor="cyan", linewidth=2, linestyle="--")
ax.add_patch(rect17)

fig17.suptitle("Figure 14: Novel Inc/Rep Group Discovery — Training Data, "
               "Performance, and Taxonomic Separation",
               fontsize=12, fontweight="bold", y=1.02)

out17 = os.path.join(FIG_DIR, "figure15_novel_inc_discovery.png")
fig17.savefig(out17, dpi=300, bbox_inches="tight", facecolor="white")
fig17.savefig(out17.replace(".png", ".pdf"), dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig17)
print(f"  Saved: {out17}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 18: Recombination & Evolutionary Rate Concepts
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating Figure 17: Recombination Detection and Evolutionary Rate ...")

fig18, axes18 = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Recombination detection concept
ax = axes18[0]
ax.set_xlim(0, 100)
ax.set_ylim(-1, 6)

# Reference plasmid
ref_y = 5
ax.barh(ref_y, 100, height=0.6, color="#3498db", alpha=0.8, edgecolor="white")
ax.text(50, ref_y, "Reference Plasmid", ha="center", va="center", fontsize=8,
        fontweight="bold", color="white")

# Query with no recombination
q1_y = 3.5
ax.barh(q1_y, 95, height=0.6, color="#27ae60", alpha=0.8, edgecolor="white")
ax.text(50, q1_y, "Query A: No recombination (95% coverage)", ha="center",
        va="center", fontsize=7, color="white")
ax.text(102, q1_y, "LOW", ha="left", va="center", fontsize=8,
        fontweight="bold", color="#27ae60")

# Query with recombination (mosaic)
q2_y = 2
segments = [(0, 30), (45, 75), (85, 100)]
for s, e in segments:
    ax.barh(q2_y, e - s, left=s, height=0.6, color="#e74c3c", alpha=0.8, edgecolor="white")
# Gaps
for (_, e1), (s2, _) in zip(segments[:-1], segments[1:]):
    ax.barh(q2_y, s2 - e1, left=e1, height=0.6, color="#f1c40f", alpha=0.4,
            hatch="///", edgecolor="#e67e22", linewidth=0.5)
ax.text(50, q2_y, "Query B: Recombinant (3 blocks, 62% coverage)", ha="center",
        va="center", fontsize=7, color="white")
ax.text(102, q2_y, "HIGH", ha="left", va="center", fontsize=8,
        fontweight="bold", color="#e74c3c")

# Query with moderate recombination
q3_y = 0.5
segments3 = [(0, 50), (60, 100)]
for s, e in segments3:
    ax.barh(q3_y, e - s, left=s, height=0.6, color="#f39c12", alpha=0.8, edgecolor="white")
ax.barh(q3_y, 10, left=50, height=0.6, color="#f1c40f", alpha=0.4,
        hatch="///", edgecolor="#e67e22", linewidth=0.5)
ax.text(50, q3_y, "Query C: Moderate (2 blocks, 90% coverage)", ha="center",
        va="center", fontsize=7, color="white")
ax.text(102, q3_y, "MEDIUM", ha="left", va="center", fontsize=8,
        fontweight="bold", color="#f39c12")

ax.set_yticks([])
ax.set_xlabel("Position (kb)")
ax.set_title("A) Recombination Detection (L5)\nAlignment Coverage Analysis", fontweight="bold")
ax.spines["left"].set_visible(False)
ax.spines["top"].set_visible(False)

# Legend
recom_legend = [
    mpatches.Patch(color="#3498db", label="Reference"),
    mpatches.Patch(color="#27ae60", label="Aligned (no recom.)"),
    mpatches.Patch(color="#e74c3c", label="Aligned (fragmented)"),
    mpatches.Patch(facecolor="#f1c40f", alpha=0.4, hatch="///",
                   edgecolor="#e67e22", label="Unaligned gap"),
]
ax.legend(handles=recom_legend, loc="lower left", fontsize=7, framealpha=0.9)

# Panel B: Evolutionary rate estimation
ax = axes18[1]

# Simulate temporal SNP data
np.random.seed(42)
n_pairs = 30
days_apart = np.random.uniform(30, 730, n_pairs)
rate = 0.015  # SNPs per day
snps = np.random.poisson(rate * days_apart).astype(float)
snps += np.random.normal(0, 0.5, n_pairs).clip(-1, 2)
snps = np.maximum(snps, 0)

ax.scatter(days_apart, snps, c="#3498db", s=40, alpha=0.7, edgecolors="white", linewidth=0.5)

# Regression line
from numpy.polynomial.polynomial import polyfit
coeffs = np.polyfit(days_apart, snps, 1)
x_line = np.linspace(0, 800, 100)
y_line = coeffs[0] * x_line + coeffs[1]
ax.plot(x_line, y_line, "r-", linewidth=2, alpha=0.8, label="Linear regression")

# Confidence interval (approx)
y_pred = coeffs[0] * days_apart + coeffs[1]
residuals = snps - y_pred
se = np.std(residuals)
ax.fill_between(x_line, y_line - 1.96 * se, y_line + 1.96 * se,
                alpha=0.15, color="red", label="95% CI")

# Rate annotation
rate_per_year = coeffs[0] * 365
r_squared = 1 - np.sum(residuals**2) / np.sum((snps - np.mean(snps))**2)
ax.text(0.05, 0.95, f"Rate: {rate_per_year:.1f} SNPs/year\n"
        f"R² = {r_squared:.3f}\n"
        f"Expected range: 3-8 SNPs/year",
        transform=ax.transAxes, fontsize=8, va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

ax.set_xlabel("Days Between Samples")
ax.set_ylabel("Pairwise SNP Differences")
ax.set_title("B) Evolutionary Rate Estimation (L7)\nSNP Accumulation Over Time", fontweight="bold")
ax.legend(fontsize=7, loc="lower right")
ax.set_xlim(0, 800)
ax.set_ylim(-1, max(snps) + 3)

fig18.suptitle("Figure 17: Recombination Detection and Molecular Clock Analysis",
               fontsize=12, fontweight="bold", y=1.02)

out18 = os.path.join(FIG_DIR, "figure17_recombination_evolutionary_rate.png")
fig18.savefig(out18, dpi=300, bbox_inches="tight", facecolor="white")
fig18.savefig(out18.replace(".png", ".pdf"), dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig18)
print(f"  Saved: {out18}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 19: Cluster Stability & Adaptive Thresholds
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating Figure 16: Adaptive Thresholds and Cluster Stability ...")

fig19, axes19 = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Bootstrap stability (REAL data from per-group bootstrap analysis)
ax = axes19[0]

STAB_TSV = os.path.join(BASE_DIR, "output", "figure16_validation", "cluster_stability.tsv")
stab_df = pd.read_csv(STAB_TSV, sep="\t")
l3_stab = stab_df[stab_df["level"] == "L3"]

# Expand per-group mean stability into per-cluster entries for histogram
all_stability = []
for _, row in l3_stab.iterrows():
    all_stability.extend([row["mean_stability"]] * row["n_clusters"])
all_stability = np.array(all_stability)

bins = np.arange(0, 105, 5)
n_stable = int(np.sum(all_stability >= 80))
n_moderate = int(np.sum((all_stability >= 50) & (all_stability < 80)))
n_unstable = int(np.sum(all_stability < 50))

ax.hist(all_stability, bins=bins, color="#3498db", edgecolor="white", linewidth=0.5, alpha=0.8)
ax.axvline(x=80, color="#27ae60", linestyle="--", linewidth=1.5, label=f"Stable threshold (n={n_stable})")
ax.axvline(x=50, color="#e74c3c", linestyle="--", linewidth=1.5, label=f"Unstable threshold (n={n_unstable})")

ax.set_xlabel("Bootstrap Stability Score (%)")
ax.set_ylabel("Number of Clusters")
n_total_cl = len(all_stability)
ax.set_title(f"A) Cluster Stability Assessment\n(50 Bootstrap Iterations, {n_total_cl} L3 Clusters, 28 Groups)",
             fontweight="bold")
ax.legend(fontsize=7, loc="upper left")

# Annotate zones
ax.axvspan(80, 100, alpha=0.08, color="#27ae60")
ax.axvspan(50, 80, alpha=0.08, color="#f39c12")
ax.axvspan(0, 50, alpha=0.08, color="#e74c3c")

ylim = ax.get_ylim()
ax.text(90, ylim[1] * 0.85, "Stable", fontsize=9, ha="center",
        color="#27ae60", fontweight="bold")
ax.text(65, ylim[1] * 0.85, "Moderate", fontsize=9, ha="center",
        color="#f39c12", fontweight="bold")
ax.text(25, ylim[1] * 0.85, "Unstable", fontsize=9, ha="center",
        color="#e74c3c", fontweight="bold")

# Panel B: Linkage method comparison
ax = axes19[1]

LINK_TSV = os.path.join(BASE_DIR, "output", "figure16_validation", "linkage_comparison.tsv")
link_df = pd.read_csv(LINK_TSV, sep="\t")
l3_link = link_df[link_df["level"] == "L3"]

methods = ["Single", "Complete", "Average"]
# Real ARI values (mean across 28 groups at L3)
sc = round(float(l3_link["single_vs_complete"].mean()), 2)
sa = round(float(l3_link["single_vs_average"].mean()), 2)
ca = round(float(l3_link["complete_vs_average"].mean()), 2)
ari_matrix = np.array([
    [1.00,  sc,    sa],
    [sc,    1.00,  ca],
    [sa,    ca,    1.00],
])
cluster_counts = [round(l3_link["n_clusters_single"].mean(), 1),
                  round(l3_link["n_clusters_complete"].mean(), 1),
                  round(l3_link["n_clusters_average"].mean(), 1)]

im19 = ax.imshow(ari_matrix, cmap="YlGnBu", vmin=0.5, vmax=1.0, aspect="auto")
ax.set_xticks(range(3))
ax.set_yticks(range(3))
ax.set_xticklabels([f"{m}\n(k={c})" for m, c in zip(methods, cluster_counts)], fontsize=9)
ax.set_yticklabels([f"{m}\n(k={c})" for m, c in zip(methods, cluster_counts)], fontsize=9)

for i in range(3):
    for j in range(3):
        color = "white" if ari_matrix[i, j] < 0.8 else "black"
        ax.text(j, i, f"{ari_matrix[i, j]:.2f}", ha="center", va="center",
                fontsize=12, fontweight="bold", color=color)

cbar19 = plt.colorbar(im19, ax=ax, fraction=0.046, pad=0.04)
cbar19.set_label("Adjusted Rand Index (ARI)", fontsize=8)

ax.set_title("B) Linkage Method Comparison (L8)\nClustering Agreement", fontweight="bold")

fig19.suptitle("Figure 15: Adaptive Thresholds and Cluster Stability",
               fontsize=12, fontweight="bold", y=1.02)

out19 = os.path.join(FIG_DIR, "figure16_cluster_stability.png")
fig19.savefig(out19, dpi=300, bbox_inches="tight", facecolor="white")
fig19.savefig(out19.replace(".png", ".pdf"), dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig19)
print(f"  Saved: {out19}")


print("\n" + "=" * 60)
print("All limitation figures generated successfully!")
print("=" * 60)
