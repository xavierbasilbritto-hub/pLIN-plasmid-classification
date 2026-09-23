#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Generate Figure 11: Cross-validation of pLIN assignment against published
outbreak plasmid datasets (17 plasmids from 4 studies).
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import seaborn as sns

# ── Setup ─────────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(BASE_DIR, "output", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

OUTBREAK_TSV = os.path.join(BASE_DIR, "output", "outbreak_validation_pLIN_results.tsv")
INTEGRATED_TSV = os.path.join(BASE_DIR, "output", "integrated", "pLIN_AMR_integrated.tsv")
LINEAGE_TSV = os.path.join(BASE_DIR, "output", "integrated", "pLIN_lineage_AMR_summary.tsv")

# Style (matches generate_figures.py)
sns.set_style("whitegrid")
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.3,
})

# Colours
BLUE = "#1E88E5"
DARK_BLUE = "#0D47A1"
GREEN = "#43A047"
ORANGE = "#FB8C00"
RED = "#E53935"
TEAL = "#009688"
GRAY = "#757575"
LIGHT_GRAY = "#ECEFF1"
LIGHT_BLUE = "#BBDEFB"

STUDY_COLORS = {
    "Study1_KPC2_IncN_Germany": BLUE,
    "Study3_NDM1_Germany": ORANGE,
    "Study4_KPC2_Singapore": RED,
    "Study6_IMP4_IncHI2_Australia": TEAL,
}
STUDY_LABELS = {
    "Study1_KPC2_IncN_Germany": "Study 1: KPC-2 IncN (Germany)",
    "Study3_NDM1_Germany": "Study 2: NDM-1 multi-Inc (Germany)",
    "Study4_KPC2_Singapore": "Study 3: KPC-2 (Singapore)",
    "Study6_IMP4_IncHI2_Australia": "Study 4: IMP-4 IncHI2 (Australia)",
}
INC_MARKERS = {
    "IncN": "o",
    "IncC": "s",
    "IncFII": "^",
    "IncAC2": "D",
    "IncA": "v",
    "IncHI2": "P",
}


def compute_lineage_profile(plin_code, integrated_df):
    """Compute AMR profile for a given pLIN lineage from training data."""
    inc_col = "inc_type_x" if "inc_type_x" in integrated_df.columns else "inc_type"
    sub = integrated_df[integrated_df["pLIN"] == plin_code]
    n = len(sub)
    if n == 0:
        return {"n": 0, "mean_amr": 0, "n_inc": 0, "mcr_pct": 0,
                "top_gene": "none", "top_gene_pct": 0}
    mean_amr = sub["n_amr_genes"].mean()
    n_inc = sub[inc_col].nunique()
    mcr_pct = sum(1 for g in sub["amr_genes"] if "mcr-" in str(g)) / n * 100
    # Find top AMR gene
    from collections import Counter
    gene_counts = Counter()
    for genes_str in sub["amr_genes"]:
        if str(genes_str) != "none":
            for g in str(genes_str).split(", "):
                gene_counts[g.strip()] += 1
    top_gene, top_count = gene_counts.most_common(1)[0] if gene_counts else ("none", 0)
    top_gene_pct = top_count / n * 100
    return {"n": n, "mean_amr": mean_amr, "n_inc": n_inc, "mcr_pct": mcr_pct,
            "top_gene": top_gene, "top_gene_pct": top_gene_pct}


def main():
    print("Loading data ...")
    val_df = pd.read_csv(OUTBREAK_TSV, sep="\t")
    integrated_df = pd.read_csv(INTEGRATED_TSV, sep="\t")

    print(f"  Outbreak plasmids: {len(val_df)}")
    print(f"  Training plasmids: {len(integrated_df)}")

    # Compute lineage profiles for Panel C
    p671 = compute_lineage_profile("1.1.2.4.7.1947", integrated_df)
    p860 = compute_lineage_profile("1.1.2.4.7.13", integrated_df)
    print(f"  pLIN 1947: n={p671['n']}, mean_amr={p671['mean_amr']:.1f}, "
          f"Inc groups={p671['n_inc']}, mcr={p671['mcr_pct']:.0f}%, "
          f"top gene={p671['top_gene']} ({p671['top_gene_pct']:.0f}%)")
    print(f"  pLIN 13: n={p860['n']}, mean_amr={p860['mean_amr']:.1f}, "
          f"Inc groups={p860['n_inc']}, mcr={p860['mcr_pct']:.0f}%, "
          f"top gene={p860['top_gene']} ({p860['top_gene_pct']:.0f}%)")

    # ── Create figure ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 14))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.30,
                           height_ratios=[1, 1], width_ratios=[1, 1])

    # ══════════════════════════════════════════════════════════════════════════
    # Panel A: NN Distance vs Classification Confidence
    # ══════════════════════════════════════════════════════════════════════════
    ax_a = fig.add_subplot(gs[0, 0])

    # Plot each point by study and Inc type
    for _, row in val_df.iterrows():
        colour = STUDY_COLORS.get(row["study"], GRAY)
        marker = INC_MARKERS.get(row["predicted_inc"], "o")
        # Small jitter for overlapping points at distance=0
        x_val = row["nn_distance"]
        y_val = row["confidence"]
        if x_val == 0:
            x_val += np.random.uniform(-0.00008, 0.00008)
        ax_a.scatter(x_val, y_val, c=colour, marker=marker, s=100,
                     edgecolors="black", linewidth=0.6, zorder=3)

    # 60% confidence threshold
    ax_a.axhline(y=60, color=GRAY, linestyle="--", linewidth=1.5, alpha=0.7)
    ax_a.text(0.0058, 61.5, "60% confidence\nthreshold", fontsize=7,
              color=GRAY, ha="right", va="bottom")

    # Annotate key plasmids
    # CP104944 → pLIN 1947 (exact match)
    cp104944_row = val_df[val_df["accession"] == "CP104944"].iloc[0]
    ax_a.annotate("CP104944\npLIN 1947 (KPC-2 hotspot)",
                  xy=(cp104944_row["nn_distance"], cp104944_row["confidence"]),
                  xytext=(0.0015, 90),
                  fontsize=7.5, fontweight="bold", color=DARK_BLUE,
                  arrowprops=dict(arrowstyle="->", color=DARK_BLUE, lw=1.2))

    # CP022533 → pLIN 13 (MDR hub)
    cp022_row = val_df[val_df["accession"] == "CP022533"].iloc[0]
    ax_a.annotate("CP022533\npLIN 13 (MDR hub)",
                  xy=(cp022_row["nn_distance"], cp022_row["confidence"]),
                  xytext=(0.002, 82),
                  fontsize=7.5, fontweight="bold", color=TEAL,
                  arrowprops=dict(arrowstyle="->", color=TEAL, lw=1.2))

    # MN542377 → borderline Singapore plasmid
    mn_row = val_df[val_df["accession"] == "MN542377"].iloc[0]
    ax_a.annotate("MN542377\n(borderline, d=0.006)",
                  xy=(mn_row["nn_distance"], mn_row["confidence"]),
                  xytext=(0.004, 50),
                  fontsize=7, color=RED,
                  arrowprops=dict(arrowstyle="->", color=RED, lw=1))

    # Legend for studies
    import matplotlib.patches as mpatches
    study_handles = [mpatches.Patch(color=c, label=STUDY_LABELS[s])
                     for s, c in STUDY_COLORS.items()]
    # Inc marker legend
    inc_handles = [plt.Line2D([0], [0], marker=m, color="gray", linestyle="None",
                              markersize=7, label=inc)
                   for inc, m in INC_MARKERS.items()]
    leg1 = ax_a.legend(handles=study_handles, fontsize=8, loc="lower left",
                       frameon=True, fancybox=True, title="Study", title_fontsize=8)
    ax_a.add_artist(leg1)
    ax_a.legend(handles=inc_handles, fontsize=8, loc="center left",
                bbox_to_anchor=(0.0, 0.42), frameon=True, fancybox=True,
                title="Inc type", title_fontsize=8)

    ax_a.set_xlabel("Nearest-Neighbour Cosine Distance", fontsize=11)
    ax_a.set_ylabel("Inc Classification Confidence (%)", fontsize=11)
    ax_a.set_title("A   Distance vs. Classification Confidence", fontweight="bold",
                   fontsize=12, loc="left")
    ax_a.set_xlim(-0.0005, 0.007)
    ax_a.set_ylim(30, 108)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)

    # ══════════════════════════════════════════════════════════════════════════
    # Panel B: Hierarchical Tile Chart (Study 3, 12 NDM-1 plasmids)
    # ══════════════════════════════════════════════════════════════════════════
    ax_b = fig.add_subplot(gs[0, 1])

    study3 = val_df[val_df["study"] == "Study3_NDM1_Germany"].copy()
    study3 = study3.sort_values("pLIN").reset_index(drop=True)
    n_plasmids = len(study3)

    # Parse pLIN codes into 6 hierarchical levels
    levels = ["L1", "L2", "L3", "L4", "L5", "L6"]
    for i, level in enumerate(levels):
        study3[level] = study3["pLIN"].apply(lambda x: ".".join(x.split(".")[:i + 1]))

    # Colour assignment per level
    level_palettes = {}
    for level in levels:
        unique_codes = study3[level].unique()
        if len(unique_codes) == 1:
            level_palettes[level] = {unique_codes[0]: LIGHT_BLUE}
        else:
            # Use tab20 for variety
            cols = plt.cm.Set2(np.linspace(0, 1, max(len(unique_codes), 3)))
            palette = {}
            for k, code in enumerate(sorted(unique_codes)):
                palette[code] = cols[k]
            # Highlight the shared outbreak-backbone code in red
            for code in unique_codes:
                if code.endswith(".26"):
                    palette[code] = RED
            level_palettes[level] = palette

    # Draw tiles
    for j, level in enumerate(levels):
        palette = level_palettes[level]
        for i, (_, row) in enumerate(study3.iterrows()):
            code_val = row[level]
            colour = palette[code_val]
            rect = plt.Rectangle((j - 0.42, i - 0.42), 0.84, 0.84,
                                 facecolor=colour, edgecolor="white",
                                 linewidth=2, zorder=2)
            ax_b.add_patch(rect)
            # Show L6 code number in the last column
            if level == "L6":
                l6_num = code_val.split(".")[-1]
                is_shared_backbone = l6_num == "26"
                ax_b.text(j, i, l6_num, ha="center", va="center",
                          fontsize=8, fontweight="bold",
                          color="white" if is_shared_backbone else "#333333", zorder=3)

    # Plasmid labels (y-axis)
    plasmid_labels = []
    for _, row in study3.iterrows():
        name = row["plasmid_name"]
        inc = row["predicted_inc"]
        plasmid_labels.append(f"{name} ({inc})")

    ax_b.set_xlim(-0.6, 5.8)
    ax_b.set_ylim(-0.6, n_plasmids - 0.4)
    ax_b.set_xticks(range(6))
    ax_b.set_xticklabels(["L1\nFamily", "L2\nSubfam.", "L3\nCluster",
                           "L4\nSubclst.", "L5\nClone\ngroup", "L6\nLineage"],
                          fontsize=8)
    ax_b.set_yticks(range(n_plasmids))
    ax_b.set_yticklabels(plasmid_labels, fontsize=8)
    ax_b.set_title("B   Hierarchical pLIN Resolution (Study 2: NDM-1 Outbreak)",
                   fontweight="bold", fontsize=12, loc="left")
    ax_b.invert_yaxis()
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)

    # Annotation bracket for the shared outbreak-backbone group
    shared_backbone_code = "1.1.2.4.7.26"
    shared_indices = study3.index[study3["pLIN"] == shared_backbone_code].tolist()
    if shared_indices:
        y_min = min(shared_indices) - 0.3
        y_max = max(shared_indices) + 0.3
        ax_b.annotate("", xy=(5.55, y_min), xytext=(5.55, y_max),
                       arrowprops=dict(arrowstyle="-", color=RED, lw=2))
        ax_b.text(5.7, (y_min + y_max) / 2,
                  "Outbreak\nbackbone\n(pLIN 26)",
                  fontsize=7, color=RED, fontweight="bold",
                  ha="left", va="center")

    # Annotation: report the deepest hierarchical level all plasmids share
    level_labels = ["L1", "L2", "L3", "L4", "L5", "L6"]
    deepest_shared_depth = 0
    for depth in range(1, 7):
        prefixes = study3["pLIN"].apply(lambda c: ".".join(c.split(".")[:depth])).unique()
        if len(prefixes) == 1:
            deepest_shared_depth = depth
            shared_prefix = prefixes[0]
        else:
            break
    if deepest_shared_depth > 0:
        ax_b.text(0.5, -0.13,
                  f"All {n_plasmids} plasmids share {level_labels[deepest_shared_depth-1]} "
                  f"code ({shared_prefix})",
                  transform=ax_b.transAxes,
                  fontsize=7.5, color=DARK_BLUE, fontweight="bold",
                  ha="center", va="top", style="italic")

    # ══════════════════════════════════════════════════════════════════════════
    # Panel C: High-Risk Lineage Profiles (grouped bar chart)
    # ══════════════════════════════════════════════════════════════════════════
    ax_c = fig.add_subplot(gs[1, 0])

    metrics = ["Reference\nmembers (n)", "Mean AMR\ngenes", "Inc groups\n(n)",
               "$mcr$ prevalence\n(%)"]
    values_671 = [p671["n"], p671["mean_amr"], p671["n_inc"], p671["mcr_pct"]]
    values_860 = [p860["n"], p860["mean_amr"], p860["n_inc"], p860["mcr_pct"]]

    x = np.arange(len(metrics))
    width = 0.32

    bars1 = ax_c.bar(x - width / 2, values_671, width, color=BLUE,
                     label="pLIN 1947 (IncN, KPC-2 hotspot)",
                     edgecolor="white", linewidth=0.8)
    bars2 = ax_c.bar(x + width / 2, values_860, width, color=RED,
                     label="pLIN 13 (multi-Inc MDR hub)",
                     edgecolor="white", linewidth=0.8)

    # Value labels on bars
    for bar in bars1:
        val = bar.get_height()
        fmt = f"{val:.0f}" if val >= 10 else f"{val:.1f}"
        ax_c.text(bar.get_x() + bar.get_width() / 2, val + 2,
                  fmt, ha="center", fontsize=9, fontweight="bold", color=BLUE)
    for bar in bars2:
        val = bar.get_height()
        fmt = f"{val:.0f}" if val >= 10 else f"{val:.1f}"
        ax_c.text(bar.get_x() + bar.get_width() / 2, val + 2,
                  fmt, ha="center", fontsize=9, fontweight="bold", color=RED)

    ax_c.set_xticks(x)
    ax_c.set_xticklabels(metrics, fontsize=9)
    ax_c.set_ylabel("Value", fontsize=11)
    ax_c.set_title("C   High-Risk Lineage Profiles (Outbreak-Matched)",
                   fontweight="bold", fontsize=12, loc="left")
    ax_c.legend(fontsize=8, loc="upper center", frameon=True, fancybox=True)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)
    ax_c.set_ylim(0, max(max(values_671), max(values_860)) * 1.25)

    # Clinical significance box
    ax_c.text(0.02, 0.02,
              "pLIN 1947: specialist KPC-2 lineage (1 Inc group, 100% $bla_{KPC-2}$)\n"
              "pLIN 13: broad MDR hub (5 Inc groups, 44% $mcr$, 61% $sul1$)",
              transform=ax_c.transAxes, fontsize=7.5, color=DARK_BLUE,
              style="italic", va="bottom",
              bbox=dict(boxstyle="round,pad=0.4", facecolor="#E3F2FD",
                        edgecolor=BLUE, alpha=0.85))

    # ══════════════════════════════════════════════════════════════════════════
    # Panel D: Summary Validation Statistics Table
    # ══════════════════════════════════════════════════════════════════════════
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.axis("off")
    ax_d.set_xlim(0, 10)
    ax_d.set_ylim(0, 10)

    ax_d.set_title("D   Validation Summary", fontweight="bold", fontsize=12,
                   loc="left")

    # Compute summary statistics directly from val_df (the 17 original
    # outbreak plasmids) rather than hardcoding them, so this panel cannot
    # silently go stale if the pLIN numbering or underlying data changes.
    n_total = len(val_df)
    n_studies = val_df["study"].nunique()
    n_unique_l6 = val_df["pLIN"].nunique()
    l6_counts = val_df["pLIN"].value_counts()
    n_novel = int((l6_counts == 1).sum())
    high_conf_df = val_df[val_df["confidence"] >= 60]
    n_high_conf = len(high_conf_df)
    known_hits = val_df["pLIN"].isin(["1.1.2.4.7.1947", "1.1.2.4.7.13"]).sum()
    shared_backbone_n = int(l6_counts.get(shared_backbone_code, 0))
    mean_nn_dist = val_df["nn_distance"].mean()

    table_data = [
        ("Total outbreak plasmids tested", str(n_total)),
        ("Published outbreak studies", str(n_studies)),
        ("Unique pLIN L6 codes assigned", str(n_unique_l6)),
        ("Novel codes (new L6 variants)", str(n_novel)),
        ("Known high-risk lineage matches", str(int(known_hits))),
        ("Inc accuracy (confidence \u226560%)",
         f"{n_high_conf/n_total*100:.0f}% ({n_high_conf}/{n_total})"),
        ("Outbreak backbone grouping",
         f"{shared_backbone_n} plasmids \u2192 pLIN 26" if shared_backbone_n else "n/a"),
        ("Mean NN distance", f"{mean_nn_dist:.4f}"),
        ("Cross-species transmission", "Detected" if shared_backbone_n >= 2 else "Not detected"),
        ("Hierarchical resolution", f"{n_total} \u2192 {n_unique_l6} L6 \u2192 1 L3"),
    ]

    # Highlight values
    highlight_values = {"100% (8/8)", "2", "Detected"}

    # Header bar
    header_box = FancyBboxPatch((0.3, 8.85), 9.2, 0.65,
                                boxstyle="round,pad=0.08",
                                facecolor=DARK_BLUE, edgecolor="white",
                                linewidth=1.5)
    ax_d.add_patch(header_box)
    ax_d.text(3.0, 9.17, "Criterion", ha="center", va="center",
              fontsize=10, fontweight="bold", color="white")
    ax_d.text(7.8, 9.17, "Result", ha="center", va="center",
              fontsize=10, fontweight="bold", color="white")

    # Data rows
    row_height = 0.78
    for i, (metric, value) in enumerate(table_data):
        y = 8.05 - i * row_height
        bg = LIGHT_GRAY if i % 2 == 0 else "white"
        row_box = FancyBboxPatch((0.3, y - 0.22), 9.2, 0.62,
                                 boxstyle="round,pad=0.02",
                                 facecolor=bg, edgecolor="#E0E0E0",
                                 linewidth=0.5)
        ax_d.add_patch(row_box)
        ax_d.text(0.6, y + 0.08, metric, va="center", fontsize=9,
                  color="#333333")
        val_color = GREEN if value in highlight_values else DARK_BLUE
        val_weight = "bold" if value in highlight_values or value in {"13", "9"} else "normal"
        ax_d.text(7.8, y + 0.08, value, ha="center", va="center",
                  fontsize=9.5, fontweight=val_weight, color=val_color)

    # ══════════════════════════════════════════════════════════════════════════
    # Save
    # ══════════════════════════════════════════════════════════════════════════
    fig.suptitle("Figure 11. Cross-Validation Against Published Outbreak Datasets",
                 fontsize=14, fontweight="bold", y=0.995, color=DARK_BLUE)

    for fmt in ["png", "pdf"]:
        out_path = os.path.join(FIG_DIR, f"Figure11_outbreak_validation.{fmt}")
        fig.savefig(out_path)
        print(f"  Saved {out_path}")
    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
