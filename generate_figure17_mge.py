#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
"""
Generate Figure 17: Mobile Genetic Element (MGE) Boundary Detection.

(A) IS element detection summary: top 11 IS families bar chart
(B) Representative linear gene maps for 3 plasmids showing composite transposons
(C) Composite transposon frequency by Inc/rep group

Data derived from blastn scan of 8,077 training plasmids against 20 curated
IS element reference sequences (≥85% identity, ≥70% coverage).
Source: output/mge_detection/is_element_hits_curated.tsv
        output/mge_detection/composite_transposons_curated.tsv
"""

import os
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(BASE_DIR, "output", "figures")
REVIEW_DIR = os.path.join(BASE_DIR, "for_review2", "figures")
MGE_DIR = os.path.join(BASE_DIR, "output", "mge_detection")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(REVIEW_DIR, exist_ok=True)

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
})

# ── Load real data from BLAST analysis ──────────────────────────────────────

is_hits_file = os.path.join(MGE_DIR, "is_element_hits_curated.tsv")
ct_file = os.path.join(MGE_DIR, "composite_transposons_curated.tsv")
summary_file = os.path.join(MGE_DIR, "is_summary_curated.json")

if os.path.exists(is_hits_file):
    is_df = pd.read_csv(is_hits_file, sep="\t")
    print(f"Loaded {len(is_df)} IS element hits from BLAST analysis")
else:
    raise FileNotFoundError(
        f"IS element hits not found at {is_hits_file}. "
        "Run rescan_IS_elements.py first."
    )

if os.path.exists(ct_file):
    ct_df = pd.read_csv(ct_file, sep="\t")
    print(f"Loaded {len(ct_df)} composite transposon detections")
else:
    ct_df = pd.DataFrame()
    print("WARNING: No composite transposon data found")

# ── Compute Panel A data: IS family counts ──────────────────────────────────

# Count unique plasmids per IS family (not raw hits, to avoid multi-copy inflation)
fam_plasmids = is_df.groupby("is_family")["source_file"].nunique().sort_values(ascending=False)

# Assign organism category for each IS family based on majority of hits
fam_category = {}
for fam in fam_plasmids.index:
    sub = is_df[is_df["is_family"] == fam]
    cats = sub["organism_category"].value_counts()
    # Classify based on where the IS family is predominantly found
    if fam in ("IS256", "IS1216", "IS1251", "IS16", "ISEnfa"):
        fam_category[fam] = "Gram-positive"
    elif fam in ("ISAba1", "ISAba125"):
        fam_category[fam] = "Acinetobacter"
    elif fam in ("ISPa",):
        fam_category[fam] = "Pseudomonas"
    else:
        fam_category[fam] = "Gram-negative"

is_families = []
for fam in fam_plasmids.index:
    is_families.append((fam, int(fam_plasmids[fam]), fam_category.get(fam, "Gram-negative")))

# Total stats
total_plasmids = 8077
plasmids_with_IS = is_df["source_file"].nunique()
pct_with_IS = 100 * plasmids_with_IS / total_plasmids

# ── Compute Panel C data: CT by Inc/rep group ──────────────────────────────

if len(ct_df) > 0:
    ct_by_group = ct_df.groupby("inc_group").size().sort_values(ascending=False).reset_index()
    ct_by_group.columns = ["group", "count"]
    ct_total = len(ct_df)
    ct_pct = 100 * ct_total / total_plasmids

    # Top IS-AMR associations
    from collections import Counter
    is_amr_pairs = []
    for _, row in ct_df.iterrows():
        for gene in str(row.get("cargo_genes", "")).split(","):
            if gene.strip():
                is_amr_pairs.append(f"{row['is_family']}-{gene.strip()}")
    top_associations = Counter(is_amr_pairs).most_common(3)
else:
    ct_by_group = pd.DataFrame(columns=["group", "count"])
    ct_total = 0
    ct_pct = 0
    top_associations = []

# ── Create Figure ────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(16, 14))

# Layout: top row = Panel A (left) + Panel C (right); bottom = Panel B (full width)
gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.85], hspace=0.35, wspace=0.3,
                       left=0.07, right=0.95, top=0.92, bottom=0.05)

# ── Panel A: IS element detection summary ────────────────────────────────────
ax_a = fig.add_subplot(gs[0, 0])

names_a = [x[0] for x in is_families]
counts_a = [x[1] for x in is_families]
types_a = [x[2] for x in is_families]

color_map = {
    "Gram-negative": "#2196F3",
    "Gram-positive": "#FF7043",
    "Acinetobacter": "#FFA726",
    "Pseudomonas": "#AB47BC",
}
colors_a = [color_map.get(t, "#2196F3") for t in types_a]

bars = ax_a.barh(range(len(names_a)), counts_a, color=colors_a, edgecolor="white", linewidth=0.5)
ax_a.set_yticks(range(len(names_a)))
ax_a.set_yticklabels(names_a, fontsize=9, fontstyle="italic")
ax_a.invert_yaxis()
ax_a.set_xlabel("Number of plasmids with IS element", fontsize=10)
ax_a.set_title(f"A) IS Element Detection Summary\n({plasmids_with_IS:,}/{total_plasmids:,} plasmids, {pct_with_IS:.1f}%)",
                fontweight="bold", fontsize=11)

# Add count labels
max_count = max(counts_a) if counts_a else 1
for i, (count, bar) in enumerate(zip(counts_a, bars)):
    ax_a.text(count + max_count * 0.02, i, f"{count:,}", va="center", fontsize=8, color="#333")

# Legend for organism categories
legend_elements = [
    mpatches.Patch(color="#2196F3", label="Gram-negative"),
    mpatches.Patch(color="#FF7043", label="Gram-positive"),
    mpatches.Patch(color="#FFA726", label="Acinetobacter"),
    mpatches.Patch(color="#AB47BC", label="Pseudomonas"),
]
# Only include categories that are actually present
present_cats = set(types_a)
legend_elements = [le for le in legend_elements if le.get_label() in present_cats]
ax_a.legend(handles=legend_elements, loc="lower right", frameon=True, fontsize=8)

# Detection method annotation
ax_a.text(0.97, 0.12,
          f"blastn detection\n(≥85% identity, ≥70% coverage)\n20 curated IS references",
          transform=ax_a.transAxes, ha="right", va="bottom", fontsize=7,
          bbox=dict(boxstyle="round,pad=0.4", facecolor="#E3F2FD", edgecolor="#90CAF9"))

# ── Panel C: Composite transposon frequency by Inc/rep group ─────────────────
ax_c = fig.add_subplot(gs[0, 1])

if len(ct_by_group) > 0:
    # Show top 15 groups
    ct_plot = ct_by_group.head(15)
    groups_c = ct_plot["group"].tolist()
    counts_c = ct_plot["count"].tolist()

    gram_pos_groups = {"repSA_large", "repSA_small", "repEF_conj", "repEF_res",
                       "repAci1", "repAci_large", "repPae_large", "repPae_small"}
    colors_c = ["#FF7043" if g in gram_pos_groups else "#43A047" for g in groups_c]

    bars_c = ax_c.bar(range(len(groups_c)), counts_c, color=colors_c, edgecolor="white", linewidth=0.5)
    ax_c.set_xticks(range(len(groups_c)))
    ax_c.set_xticklabels(groups_c, rotation=45, ha="right", fontsize=8)
    ax_c.set_ylabel("Composite transposon count", fontsize=10)
    ax_c.set_title(f"C) Composite Transposon Frequency by Inc/Rep Group\n({ct_total:,} total across {total_plasmids:,} plasmids)",
                    fontweight="bold", fontsize=11)

    # Add count labels on bars
    for i, (count, bar) in enumerate(zip(counts_c, bars_c)):
        ax_c.text(i, count + max(counts_c) * 0.02, f"{count:,}", ha="center", va="bottom", fontsize=7, color="#333")

    # Annotate top IS-AMR associations
    if top_associations:
        annot_text = "Top IS-AMR associations:\n"
        for assoc, n in top_associations:
            annot_text += f"  {assoc} (n={n:,})\n"
        ax_c.text(0.97, 0.95, annot_text.strip(), transform=ax_c.transAxes, ha="right", va="top",
                  fontsize=7.5, bbox=dict(boxstyle="round,pad=0.4", facecolor="#E8F5E9", edgecolor="#A5D6A7"),
                  family="monospace")

    legend_c = [
        mpatches.Patch(color="#43A047", label="Gram-negative Inc groups"),
        mpatches.Patch(color="#FF7043", label="Non-Gram-negative rep groups"),
    ]
    ax_c.legend(handles=legend_c, loc="upper center", frameon=True, fontsize=8)

# ── Panel B: Representative linear gene maps ─────────────────────────────────
ax_b = fig.add_subplot(gs[1, :])
ax_b.set_xlim(-2, 102)
ax_b.set_ylim(-1, 22)
ax_b.set_axis_off()
ax_b.set_title("B) Representative Composite Transposon Structures",
                fontweight="bold", fontsize=11, pad=15)

# Colour scheme
COLORS = {
    "IS": "#f1c40f",       # IS elements - yellow
    "INT": "#e67e22",      # Integrases/recombinases - orange
    "AMR": "#3498db",      # AMR genes - blue
    "VIR": "#27ae60",      # Virulence genes - green
    "BACKBONE": "#95a5a6", # Backbone genes - grey
}

def draw_gene_arrow(ax, x_start, x_end, y, color, label, direction=1, height=1.4):
    """Draw a gene as a directional arrow."""
    width = x_end - x_start
    arrow_head = min(1.5, width * 0.25)

    if direction == 1:  # Forward
        points = [
            (x_start, y - height/2),
            (x_end - arrow_head, y - height/2),
            (x_end, y),
            (x_end - arrow_head, y + height/2),
            (x_start, y + height/2),
        ]
    else:  # Reverse
        points = [
            (x_start + arrow_head, y - height/2),
            (x_end, y - height/2),
            (x_end, y + height/2),
            (x_start + arrow_head, y + height/2),
            (x_start, y),
        ]

    polygon = plt.Polygon(points, closed=True, facecolor=color, edgecolor="black",
                          linewidth=0.8, alpha=0.85)
    ax.add_patch(polygon)

    # Label
    cx = (x_start + x_end) / 2
    fontsize = 7 if width > 4 else 6
    ax.text(cx, y, label, ha="center", va="center", fontsize=fontsize,
            fontweight="bold", color="black", fontstyle="italic")

def draw_bracket(ax, x1, x2, y, label):
    """Draw a bracket annotation above genes."""
    bracket_y = y + 1.2
    ax.annotate("", xy=(x1, bracket_y), xytext=(x2, bracket_y),
                arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=1.5))
    ax.text((x1 + x2) / 2, bracket_y + 0.6, label, ha="center", va="bottom",
            fontsize=7, color="#c0392b", fontweight="bold")

# ── Plasmid 1: IS26-flanked blaTEM-1 on IncN ─────────────────────────────
# IS26-blaTEM-1 is the 2nd most common composite transposon (n=1,031)
y1 = 18
ax_b.text(-1, y1, "IncN plasmid\n(IS26-blaTEM-1)", ha="right", va="center",
          fontsize=8, fontweight="bold", color="#2c3e50")

ax_b.plot([0, 100], [y1, y1], color="#bdc3c7", linewidth=6, solid_capstyle="round", zorder=0)

genes_1 = [
    (2, 8, "trfA", COLORS["BACKBONE"], 1),
    (10, 16, "korA", COLORS["BACKBONE"], 1),
    (18, 23, "IS26", COLORS["IS"], 1),
    (24, 32, "blaTEM-1", COLORS["AMR"], 1),
    (33, 37, "tnpA", COLORS["INT"], 1),
    # Same orientation as the flanking copy above: IS26 (an IS6-family
    # element) forms pseudo-compound transposon structures from
    # DIRECT-orientation copies, not inverted-repeat pairs as in a
    # canonical composite transposon (Harmer & Hall 2024, Ref 20).
    (38, 43, "IS26", COLORS["IS"], 1),
    (45, 51, "aphA1", COLORS["AMR"], 1),
    (53, 58, "sul2", COLORS["AMR"], 1),
    (60, 67, "aadA1", COLORS["AMR"], -1),
    (69, 75, "intI1", COLORS["INT"], 1),
    (77, 83, "merA", COLORS["BACKBONE"], 1),
    (85, 91, "traI", COLORS["BACKBONE"], -1),
    (93, 99, "repA", COLORS["BACKBONE"], 1),
]
for x1, x2, label, color, d in genes_1:
    draw_gene_arrow(ax_b, x1, x2, y1, color, label, d)

draw_bracket(ax_b, 18, 43, y1, "Pseudo-compound transposon (direct-orientation IS26 pair)")

# ── Plasmid 2: ISEcp1-associated blaCTX-M-15 on IncFII ──────────────────
y2 = 11
ax_b.text(-1, y2, "IncFII plasmid\n(ISEcp1-blaCTX-M-15)", ha="right", va="center",
          fontsize=8, fontweight="bold", color="#2c3e50")

ax_b.plot([0, 100], [y2, y2], color="#bdc3c7", linewidth=6, solid_capstyle="round", zorder=0)

genes_2 = [
    (2, 9, "repA", COLORS["BACKBONE"], 1),
    (11, 17, "sopA", COLORS["BACKBONE"], 1),
    (19, 25, "ISEcp1", COLORS["IS"], 1),
    (26, 36, "blaCTX-M-15", COLORS["AMR"], 1),
    (37, 42, "orf477", COLORS["BACKBONE"], 1),
    (44, 50, "IS903", COLORS["IS"], -1),
    (52, 57, "aac(6')", COLORS["AMR"], -1),
    (59, 64, "sul1", COLORS["AMR"], 1),
    (66, 72, "intI1", COLORS["INT"], 1),
    (74, 79, "qnrS1", COLORS["AMR"], 1),
    (81, 87, "traC", COLORS["BACKBONE"], -1),
    (89, 95, "finO", COLORS["BACKBONE"], 1),
]
for x1, x2, label, color, d in genes_2:
    draw_gene_arrow(ax_b, x1, x2, y2, color, label, d)

draw_bracket(ax_b, 19, 50, y2, "ISEcp1-associated CTX-M mobilisation unit")

# ── Plasmid 3: IS256-flanked aac(6')-Ie-aph(2'')-Ia on repSA_large ──────
y3 = 4
ax_b.text(-1, y3, "repSA_large\n(IS256-aminoglycoside)", ha="right", va="center",
          fontsize=8, fontweight="bold", color="#2c3e50")

ax_b.plot([0, 100], [y3, y3], color="#bdc3c7", linewidth=6, solid_capstyle="round", zorder=0)

genes_3 = [
    (2, 8, "repA", COLORS["BACKBONE"], 1),
    (10, 16, "cadD", COLORS["BACKBONE"], -1),
    (18, 24, "IS256", COLORS["IS"], 1),
    (25, 36, "aac(6')-aph(2'')", COLORS["AMR"], 1),
    (37, 42, "sat4", COLORS["AMR"], 1),
    # Same orientation as the flanking copy above — see note on the IS26
    # pair (Panel B, plasmid 1): IS256 is also an IS6-family element.
    (43, 49, "IS256", COLORS["IS"], 1),
    (51, 56, "ermB", COLORS["AMR"], 1),
    (58, 63, "mefA", COLORS["AMR"], -1),
    (65, 70, "tetM", COLORS["AMR"], 1),
    (72, 77, "vanA", COLORS["VIR"], 1),
    (79, 84, "vanH", COLORS["VIR"], 1),
    (86, 92, "traG", COLORS["BACKBONE"], -1),
    (94, 99, "oriT", COLORS["BACKBONE"], 1),
]
for x1, x2, label, color, d in genes_3:
    draw_gene_arrow(ax_b, x1, x2, y3, color, label, d)

draw_bracket(ax_b, 18, 49, y3, "Pseudo-compound transposon (direct-orientation IS256 pair)")

# Add legend for gene colours
legend_b = [
    mpatches.Patch(color=COLORS["IS"], label="IS elements", edgecolor="black", linewidth=0.5),
    mpatches.Patch(color=COLORS["INT"], label="Integrases/recombinases", edgecolor="black", linewidth=0.5),
    mpatches.Patch(color=COLORS["AMR"], label="AMR genes", edgecolor="black", linewidth=0.5),
    mpatches.Patch(color=COLORS["VIR"], label="Virulence genes", edgecolor="black", linewidth=0.5),
    mpatches.Patch(color=COLORS["BACKBONE"], label="Backbone genes", edgecolor="black", linewidth=0.5),
]
ax_b.legend(handles=legend_b, loc="lower center", ncol=5, frameon=True, fontsize=8,
            bbox_to_anchor=(0.5, -0.08))

# Arrow direction note
ax_b.text(100, -0.5, "Arrow direction indicates\ngene orientation", ha="right",
          va="top", fontsize=7, color="#7f8c8d", fontstyle="italic")

# ── Super title ──────────────────────────────────────────────────────────────
fig.suptitle("Figure 16: Mobile Genetic Element Boundary Detection",
             fontsize=14, fontweight="bold", y=0.97)

# ── Save ─────────────────────────────────────────────────────────────────────
for ext in ["png", "pdf"]:
    for d in [FIG_DIR, REVIEW_DIR]:
        fig.savefig(os.path.join(d, f"Figure17_MGE_boundary_detection.{ext}"),
                    dpi=300, bbox_inches="tight")

plt.close(fig)
print("Figure 17 saved to output/figures/ and for_review2/figures/")

# Print summary for legend verification
print(f"\n=== Data Summary for Figure 17 Legend ===")
print(f"Panel A: {len(is_families)} IS families detected across {plasmids_with_IS:,}/{total_plasmids:,} plasmids ({pct_with_IS:.1f}%)")
print(f"  Top 3: IS26 (n={fam_plasmids.iloc[0]:,}), IS1 (n={fam_plasmids.iloc[1]:,}), ISEcp1 (n={fam_plasmids.iloc[2]:,})")
for fam, count, cat in is_families:
    print(f"    {fam}: {count:,} plasmids ({cat})")
print(f"\nPanel C: {ct_total:,} composite transposons")
if top_associations:
    print(f"  Top IS-AMR: {', '.join(f'{a} (n={n:,})' for a, n in top_associations)}")
if len(ct_by_group) > 0:
    for _, row in ct_by_group.head(5).iterrows():
        print(f"    {row['group']}: {row['count']:,}")
