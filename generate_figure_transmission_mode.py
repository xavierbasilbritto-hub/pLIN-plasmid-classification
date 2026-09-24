#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
"""
Generate Figure: Combined chromosomal-plasmid typing — transmission mode
discrimination in published outbreaks.

Panel A: Bipartite network (MLST ST ↔ pLIN L6 code)
Panel B: Stacked bar chart by study (transmission mode distribution)
Panel C: Summary validation metrics
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from collections import Counter, defaultdict

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)

FIG_DIR = os.path.join(BASE_DIR, "output", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ── Color scheme ─────────────────────────────────────────────────────────────
MODE_COLORS = {
    "Clonal spread": "#E53935",
    "Horizontal plasmid transfer": "#FB8C00",
    "Same strain, different plasmids": "#1E88E5",
    "Independent": "#78909C",
}

SPECIES_COLORS = {
    "Escherichia coli": "#1565C0",
    "Klebsiella pneumoniae": "#C62828",
    "Enterobacter cloacae": "#2E7D32",
    "Citrobacter freundii": "#6A1B9A",
    "Enterobacter hormaechei": "#00838F",
    "Salmonella enterica": "#EF6C00",
    "Acinetobacter baumannii": "#4E342E",
}


def load_validation_data():
    """Load retrospective validation results."""
    study_path = os.path.join(BASE_DIR, "output",
                               "retrospective_validation_study_results.tsv")
    enriched_path = os.path.join(BASE_DIR, "output",
                                  "retrospective_host_plasmid_validation.tsv")
    study_df = pd.read_csv(study_path, sep="\t")
    enriched_df = pd.read_csv(enriched_path, sep="\t")
    return study_df, enriched_df


def shorten_study(name):
    """Shorten study names for display."""
    parts = name.split("_")
    if len(parts) >= 3:
        author = parts[0]
        year = parts[1]
        gene = parts[2] if len(parts) > 2 else ""
        return f"{author} {year}\n{gene}"
    return name


def plot_panel_a(ax, enriched_df, study_df):
    """Panel A: Bipartite network showing ST ↔ pLIN connections."""
    ax.set_title("A. MLST ST — pLIN L6 Bipartite Network", fontsize=12,
                 fontweight="bold", pad=10)

    # Filter to studies with multi-ST data
    multi_studies = study_df[study_df["n_unique_STs"] > 1]["study"].tolist()
    sub = enriched_df[enriched_df["study"].isin(multi_studies)].copy()
    sub = sub[sub["MLST_ST"].notna() & (sub["MLST_ST"] != "")]

    if len(sub) == 0:
        ax.text(0.5, 0.5, "No multi-ST data", ha="center", va="center")
        ax.axis("off")
        return

    # Get unique STs and pLINs
    sts = sorted(sub["MLST_ST"].unique())
    plins = sorted(sub["pLIN"].unique())

    # Position ST nodes on left, pLIN nodes on right
    n_st = len(sts)
    n_plin = len(plins)
    st_y = {st: i for i, st in enumerate(sts)}
    plin_y = {p: i * (n_st / max(n_plin, 1)) for i, p in enumerate(plins)}

    # Draw edges
    for _, row in sub.iterrows():
        st = row["MLST_ST"]
        plin = row["pLIN"]
        species = row.get("host_species", "")
        color = SPECIES_COLORS.get(species, "#999999")
        ax.plot([0.2, 0.8], [st_y[st], plin_y[plin]],
                color=color, alpha=0.4, linewidth=1.5)

    # Draw ST nodes (left)
    for st, y in st_y.items():
        st_species = sub[sub["MLST_ST"] == st]["host_species"].mode()
        color = SPECIES_COLORS.get(st_species.iloc[0] if len(st_species) > 0 else "", "#999")
        ax.scatter(0.2, y, s=120, c=color, edgecolors="black", linewidth=0.5, zorder=5)
        ax.text(0.15, y, st, ha="right", va="center", fontsize=7, fontweight="bold")

    # Draw pLIN nodes (right)
    for plin, y in plin_y.items():
        ax.scatter(0.8, y, s=80, c="#43A047", edgecolors="black",
                   linewidth=0.5, zorder=5, marker="s")
        plin_short = str(plin).split(".")[-1] if "." in str(plin) else str(plin)
        ax.text(0.85, y, f"pLIN {plin_short}", ha="left", va="center", fontsize=6)

    # Labels
    ax.text(0.2, n_st + 0.5, "MLST ST", ha="center", fontsize=9, fontweight="bold")
    ax.text(0.8, n_st + 0.5, "pLIN L6", ha="center", fontsize=9, fontweight="bold")

    ax.set_xlim(0, 1)
    ax.set_ylim(-1, n_st + 1.5)
    ax.axis("off")

    # Species legend
    species_in_data = sorted(set(sub["host_species"].dropna()))
    legend_handles = [mpatches.Patch(color=SPECIES_COLORS.get(sp, "#999"), label=sp)
                      for sp in species_in_data[:5]]
    ax.legend(handles=legend_handles, loc="lower left", fontsize=6,
              framealpha=0.8, title="Species", title_fontsize=7)


def plot_panel_b(ax, study_df):
    """Panel B: Stacked bar chart of transmission modes per study."""
    ax.set_title("B. Transmission Mode by Study", fontsize=12,
                 fontweight="bold", pad=10)

    studies = study_df[study_df["total_pairs"] > 0].sort_values(
        "total_pairs", ascending=True).copy()

    if len(studies) == 0:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        return

    labels = [shorten_study(s) for s in studies["study"]]
    clonal = studies["clonal_pairs"].values
    hgt = studies["hgt_pairs"].values
    other = studies["total_pairs"].values - clonal - hgt

    y_pos = np.arange(len(labels))
    bar_height = 0.6

    ax.barh(y_pos, clonal, bar_height, label="Clonal spread",
            color=MODE_COLORS["Clonal spread"], edgecolor="white", linewidth=0.5)
    ax.barh(y_pos, hgt, bar_height, left=clonal,
            label="HGT", color=MODE_COLORS["Horizontal plasmid transfer"],
            edgecolor="white", linewidth=0.5)
    ax.barh(y_pos, other, bar_height, left=clonal + hgt,
            label="Other", color=MODE_COLORS["Independent"],
            edgecolor="white", linewidth=0.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Number of Pairs", fontsize=9)
    ax.legend(fontsize=7, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_panel_c(ax, study_df):
    """Panel C: Summary validation metrics table."""
    ax.set_title("C. Validation Summary", fontsize=12,
                 fontweight="bold", pad=10)
    ax.axis("off")

    n_studies = len(study_df)
    concordant = study_df["concordant"].sum()
    multi_st_correct = study_df["multi_ST_correct"].sum()
    clonal_correct = study_df["clonal_correct"].sum()
    hgt_detected = (study_df["hgt_pairs"] > 0).sum()
    clonal_detected = (study_df["clonal_pairs"] > 0).sum()
    total_species = study_df["n_species"].max()

    metrics = [
        ("Studies evaluated", f"{n_studies}"),
        ("Host species", f"7"),
        ("Unique MLST STs", f"20"),
        ("Overall concordance", f"{concordant}/{n_studies} ({100*concordant/n_studies:.1f}%)"),
        ("Multi-ST detection", f"{multi_st_correct}/{n_studies} ({100*multi_st_correct/n_studies:.1f}%)"),
        ("Clonal detection", f"{clonal_correct}/{n_studies} ({100*clonal_correct/n_studies:.1f}%)"),
        ("Studies with HGT pairs", f"{hgt_detected}"),
        ("Studies with clonal pairs", f"{clonal_detected}"),
    ]

    table_data = [[m[0], m[1]] for m in metrics]
    table = ax.table(cellText=table_data, colLabels=["Metric", "Value"],
                     loc="center", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.4)

    # Style header
    for j in range(2):
        table[(0, j)].set_facecolor("#1565C0")
        table[(0, j)].set_text_props(color="white", fontweight="bold")

    # Alternate row colors
    for i in range(1, len(metrics) + 1):
        for j in range(2):
            if i % 2 == 0:
                table[(i, j)].set_facecolor("#E3F2FD")
            else:
                table[(i, j)].set_facecolor("white")


def main():
    print("Generating transmission mode figure...")
    study_df, enriched_df = load_validation_data()

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3,
                          left=0.08, right=0.95, top=0.92, bottom=0.05)

    # Panel A: Bipartite network (top-left, spans 2 rows)
    ax_a = fig.add_subplot(gs[:, 0])
    plot_panel_a(ax_a, enriched_df, study_df)

    # Panel B: Stacked bar (top-right)
    ax_b = fig.add_subplot(gs[0, 1])
    plot_panel_b(ax_b, study_df)

    # Panel C: Summary table (bottom-right)
    ax_c = fig.add_subplot(gs[1, 1])
    plot_panel_c(ax_c, study_df)

    fig.suptitle("Figure 12. Combined Chromosomal-Plasmid Typing: Transmission Mode Discrimination",
                 fontsize=14, fontweight="bold", y=0.97)

    # Save
    png_path = os.path.join(FIG_DIR, "Figure13_transmission_mode.png")
    pdf_path = os.path.join(FIG_DIR, "Figure13_transmission_mode.pdf")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


if __name__ == "__main__":
    main()
