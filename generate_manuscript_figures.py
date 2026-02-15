#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Generate methods and architecture figures for the pLIN manuscript.
Produces Figures 8–12: pipeline overview, classification system, comparison
table, v2.1 upgrade overview, and GUI screenshot diagram.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe

# ── Setup ─────────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(BASE_DIR, "output", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.2,
})

# Colors
BLUE = "#1E88E5"
DARK_BLUE = "#0D47A1"
GREEN = "#43A047"
ORANGE = "#FB8C00"
RED = "#E53935"
PURPLE = "#8E24AA"
TEAL = "#009688"
GRAY = "#757575"
LIGHT_GRAY = "#ECEFF1"
LIGHT_BLUE = "#BBDEFB"


def add_fancy_box(ax, x, y, w, h, color, text, fontsize=9, text_color="white",
                  alpha=1.0, bold=True, rounded=True):
    """Add a rounded rectangle with centered text."""
    style = "round,pad=0.1" if rounded else "square,pad=0.05"
    box = FancyBboxPatch((x, y), w, h, boxstyle=style,
                         facecolor=color, edgecolor="white", linewidth=1.5, alpha=alpha)
    ax.add_patch(box)
    weight = "bold" if bold else "normal"
    ax.text(x + w/2, y + h/2, text, ha="center", va="center",
            fontsize=fontsize, color=text_color, fontweight=weight,
            path_effects=[pe.withStroke(linewidth=0.5, foreground="black")] if text_color == "white" else [])


def add_arrow(ax, x1, y1, x2, y2, color="#555555"):
    """Add an arrow between two points."""
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=2,
                                connectionstyle="arc3,rad=0"))


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 8: pLIN Pipeline Overview (Methods Figure)
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 8: Pipeline Overview ...")

fig8, ax = plt.subplots(1, 1, figsize=(16, 9))
ax.set_xlim(0, 16)
ax.set_ylim(0, 9)
ax.axis("off")
ax.set_aspect("equal")

# Title
ax.text(8, 8.6, "Figure 8. pLIN Analysis Pipeline Overview",
        ha="center", fontsize=14, fontweight="bold", color=DARK_BLUE)
ax.text(8, 8.25, "End-to-end workflow from FASTA input to classified plasmids with AMR surveillance",
        ha="center", fontsize=10, color=GRAY)

# ── Row 1: Input ──
add_fancy_box(ax, 0.3, 6.8, 2.5, 1.0, BLUE, "INPUT\n\nPlasmid FASTA\n(.fasta/.fa/.fna)", fontsize=9)
add_arrow(ax, 2.9, 7.3, 3.4, 7.3)

# Step 1: 4-mer computation
add_fancy_box(ax, 3.5, 6.8, 2.8, 1.0, ORANGE, "STEP 1\n\n4-mer Frequency\nVectors (256D)", fontsize=9)
add_arrow(ax, 6.4, 7.3, 6.9, 7.3)

# Step 2: Distance matrix
add_fancy_box(ax, 7.0, 6.8, 2.8, 1.0, GREEN, "STEP 2\n\nPairwise Cosine\nDistance Matrix", fontsize=9)
add_arrow(ax, 9.9, 7.3, 10.4, 7.3)

# Step 3: Clustering
add_fancy_box(ax, 10.5, 6.8, 2.8, 1.0, PURPLE, "STEP 3\n\nSingle-Linkage\nClustering", fontsize=9)
add_arrow(ax, 13.4, 7.3, 13.9, 7.3)

# Step 4: pLIN codes
add_fancy_box(ax, 14.0, 6.8, 1.7, 1.0, DARK_BLUE, "STEP 4\n\npLIN\nCodes", fontsize=9)

# ── Row 2: Inc detection + AMR (parallel paths) ──
add_arrow(ax, 1.55, 6.8, 1.55, 6.3)

# Inc Group Detection (branch from input)
add_fancy_box(ax, 0.3, 5.2, 2.5, 1.0, "#FF7043", "STEP 5\n\nKNN Inc Group\nDetection (k=5)", fontsize=9)
add_arrow(ax, 2.9, 5.7, 3.5, 5.7)

# AMRFinderPlus
add_fancy_box(ax, 3.5, 5.2, 2.8, 1.0, RED, "STEP 6\n\nAMRFinderPlus\nAMR/Stress/Virulence", fontsize=9)
add_arrow(ax, 6.4, 5.7, 6.9, 5.7)

# MOBsuite
add_fancy_box(ax, 7.0, 5.2, 2.8, 1.0, TEAL, "STEP 7\n\nMOBsuite\nMobility Typing", fontsize=9)
add_arrow(ax, 9.9, 5.7, 10.4, 5.7)

# Integration
add_fancy_box(ax, 10.5, 5.2, 2.8, 1.0, "#5C6BC0", "STEP 8\n\nIntegration\npLIN + AMR + Mobility", fontsize=9)
add_arrow(ax, 13.4, 5.7, 13.9, 5.7)

# Outbreak detection
add_fancy_box(ax, 14.0, 5.2, 1.7, 1.0, RED, "STEP 8b\n\nOutbreak\nDetection", fontsize=8)

# ── Row 3: Phase 2 enhancements ──
# Mash
add_fancy_box(ax, 0.3, 3.5, 2.5, 1.0, BLUE, "STEP 9\n\nMash ANI\n(MinHash k=21)", fontsize=9)
add_arrow(ax, 2.9, 4.0, 3.5, 4.0)

# FastANI
add_fancy_box(ax, 3.5, 3.5, 2.8, 1.0, GREEN, "STEP 10\n\nFastANI\nTrue ANI", fontsize=9)
add_arrow(ax, 6.4, 4.0, 6.9, 4.0)

# SNP sub-typing
add_fancy_box(ax, 7.0, 3.5, 2.8, 1.0, ORANGE, "STEP 11\n\nminimap2 SNP\nSub-typing (L6)", fontsize=9)
add_arrow(ax, 9.9, 4.0, 10.4, 4.0)

# Temporal
add_fancy_box(ax, 10.5, 3.5, 2.8, 1.0, RED, "STEP 12\n\nTemporal Outbreak\nClustering (30d)", fontsize=9)
add_arrow(ax, 13.4, 4.0, 13.9, 4.0)

# CRISPR
add_fancy_box(ax, 14.0, 3.5, 1.7, 1.0, "#00695C", "CRISPR\nHost\nInference", fontsize=8)

# ── Row 4: Outputs ──
add_fancy_box(ax, 0.3, 1.8, 3.5, 1.0, LIGHT_GRAY, "OUTPUTS\n\nTSV Tables | PNG/PDF Figures\nZIP Bundle | JSON Reports",
              fontsize=9, text_color=DARK_BLUE)

# GUI tabs
tabs = ["Overview", "Results", "Cladogram", "AMR", "Epi", "CRISPR", "Buddy", "Export"]
tab_colors = [BLUE, GREEN, ORANGE, RED, TEAL, "#00695C", PURPLE, GRAY]
for i, (tab, tc) in enumerate(zip(tabs, tab_colors)):
    add_fancy_box(ax, 4.2 + i * 1.5, 1.8, 1.35, 1.0, tc, tab, fontsize=7)

# Phase labels
ax.text(8, 6.2, "Core Pipeline (Steps 1-8)", ha="center", fontsize=8,
        color=GRAY, style="italic")
ax.text(8, 3.2, "Phase 2 Enhancements (Steps 9-12)", ha="center", fontsize=8,
        color=GRAY, style="italic")
ax.text(8, 1.5, "Streamlit Web GUI — 8 Interactive Tabs", ha="center", fontsize=9,
        color=DARK_BLUE, fontweight="bold")

# Vertical arrows connecting rows
add_arrow(ax, 14.85, 6.8, 14.85, 6.3)
add_arrow(ax, 1.55, 5.2, 1.55, 4.6)

for fmt in ["png", "pdf"]:
    fig8.savefig(os.path.join(FIG_DIR, f"Figure8_pipeline_overview.{fmt}"))
plt.close(fig8)
print("  Saved Figure8_pipeline_overview.png/pdf")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 9: pLIN Hierarchical Classification System
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 9: Classification System ...")

fig9, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(14, 7), gridspec_kw={"width_ratios": [1, 1.2]})

# Panel A: Threshold table as visual
ax = ax_left
ax.axis("off")
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)

ax.text(5, 9.5, "A. pLIN Hierarchical Levels", ha="center", fontsize=13,
        fontweight="bold", color=DARK_BLUE)

levels = [
    ("L1", "A", "Family", 0.150, "~85%", RED),
    ("L2", "B", "Subfamily", 0.100, "~90%", ORANGE),
    ("L3", "C", "Cluster", 0.050, "~95%", "#FDD835"),
    ("L4", "D", "Subcluster", 0.020, "~98%", GREEN),
    ("L5", "E", "Clone", 0.010, "~99%", BLUE),
    ("L6", "F", "Strain", 0.001, "~99.9%", PURPLE),
]

# Column headers
headers = ["Level", "Bin", "Name", "d ≤", "ANI"]
x_positions = [0.8, 2.2, 3.5, 5.5, 7.5]
for hdr, xp in zip(headers, x_positions):
    ax.text(xp, 8.8, hdr, ha="center", fontsize=10, fontweight="bold", color="white",
            bbox=dict(boxstyle="round,pad=0.3", facecolor=DARK_BLUE, edgecolor="none"))

for i, (level, bn, name, thresh, ani, color) in enumerate(levels):
    y = 7.8 - i * 1.1
    bg = LIGHT_GRAY if i % 2 == 0 else "white"
    box = FancyBboxPatch((0.2, y - 0.35), 9.0, 0.8, boxstyle="round,pad=0.05",
                         facecolor=bg, edgecolor="#E0E0E0", linewidth=0.5)
    ax.add_patch(box)

    # Color indicator
    circle = plt.Circle((0.8, y + 0.05), 0.18, color=color, zorder=3)
    ax.add_patch(circle)
    ax.text(0.8, y + 0.05, level, ha="center", va="center", fontsize=7,
            fontweight="bold", color="white", zorder=4)

    ax.text(2.2, y + 0.05, bn, ha="center", va="center", fontsize=10, fontweight="bold", color=color)
    ax.text(3.5, y + 0.05, name, ha="center", va="center", fontsize=10, color="#333333")
    ax.text(5.5, y + 0.05, f"≤ {thresh:.3f}", ha="center", va="center", fontsize=10,
            fontfamily="monospace", color="#333333")
    ax.text(7.5, y + 0.05, ani, ha="center", va="center", fontsize=10, color="#333333")

# Example pLIN code at bottom
box_ex = FancyBboxPatch((0.5, 0.5), 8.5, 1.3, boxstyle="round,pad=0.15",
                        facecolor="#E3F2FD", edgecolor=BLUE, linewidth=2)
ax.add_patch(box_ex)
ax.text(4.75, 1.5, "Example pLIN Code", ha="center", fontsize=10, fontweight="bold", color=DARK_BLUE)
ax.text(4.75, 0.95, "1 . 1 . 3 . 5 . 12 . 45", ha="center", fontsize=18,
        fontweight="bold", fontfamily="monospace", color=DARK_BLUE)

# Panel B: Algorithm flowchart
ax = ax_right
ax.axis("off")
ax.set_xlim(0, 12)
ax.set_ylim(0, 10)

ax.text(6, 9.5, "B. Classification Algorithm", ha="center", fontsize=13,
        fontweight="bold", color=DARK_BLUE)

steps = [
    ("1", "Parse FASTA Sequences", "BioPython SeqIO reader", GREEN),
    ("2", "Compute 4-mer Vectors", "256-dim normalised frequency", ORANGE),
    ("3", "Cosine Distance Matrix", "Pairwise d(i,j) = 1 - cos(Vi,Vj)", BLUE),
    ("4", "Single-Linkage Clustering", "scipy.cluster.hierarchy.linkage", PURPLE),
    ("5", "Cut at 6 Thresholds", "fcluster at L1–L6 distances", RED),
    ("6", "Assign pLIN Codes", "Concatenate: A.B.C.D.E.F", DARK_BLUE),
]

for i, (num, title, desc, color) in enumerate(steps):
    y = 8.5 - i * 1.35
    # Number circle
    circle = plt.Circle((1.5, y), 0.35, color=color, zorder=3)
    ax.add_patch(circle)
    ax.text(1.5, y, num, ha="center", va="center", fontsize=12,
            fontweight="bold", color="white", zorder=4)
    # Title and description
    ax.text(2.3, y + 0.15, title, va="center", fontsize=11, fontweight="bold", color=color)
    ax.text(2.3, y - 0.25, desc, va="center", fontsize=9, color=GRAY)
    # Connecting line
    if i < len(steps) - 1:
        ax.plot([1.5, 1.5], [y - 0.35, y - 1.0], color="#BBDEFB", lw=2, zorder=1)

# Key properties box
props_y = 0.3
box_props = FancyBboxPatch((0.5, props_y), 11.0, 1.0, boxstyle="round,pad=0.1",
                           facecolor="#E8F5E9", edgecolor=GREEN, linewidth=1.5)
ax.add_patch(box_props)
ax.text(6, props_y + 0.7, "Key Properties", fontsize=10, fontweight="bold",
        ha="center", color=GREEN)
ax.text(6, props_y + 0.3, "Hierarchical   |   Permanent   |   Reference-Free   |   Composition-Based   |   Scalable",
        fontsize=9, ha="center", color="#333333")

fig9.tight_layout()
for fmt in ["png", "pdf"]:
    fig9.savefig(os.path.join(FIG_DIR, f"Figure9_classification_system.{fmt}"))
plt.close(fig9)
print("  Saved Figure9_classification_system.png/pdf")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 10: Comparison with Existing Methods (Table Figure)
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 10: Method Comparison ...")

fig10, ax = plt.subplots(1, 1, figsize=(14, 6))
ax.axis("off")
ax.set_xlim(0, 14)
ax.set_ylim(0, 8)

ax.text(7, 7.5, "Figure 10. Comparison of Plasmid Classification Systems",
        ha="center", fontsize=13, fontweight="bold", color=DARK_BLUE)

# Table data
methods = ["PlasmidFinder", "pMLST", "MOB-suite", "COPLA/PTU", "mge-cluster", "pLIN"]
criteria = ["Hierarchical", "Permanent\nCodes", "Reference-\nFree", "Broad\nScope", "AMR\nIntegration",
            "ANI\nValidation", "GUI"]

# Score matrix (1=yes, 0.5=partial, 0=no)
scores = [
    [0,   1,   0,   0,   0,   0,   0],    # PlasmidFinder
    [0,   1,   0,   0,   0,   0,   0],    # pMLST
    [0,   0,   0,   0.5, 0,   0,   0],    # MOB-suite
    [0.5, 0,   0,   0.5, 0,   0,   0],    # COPLA
    [0,   0,   1,   1,   0,   0,   0],    # mge-cluster
    [1,   1,   1,   1,   1,   1,   1],    # pLIN
]

# Draw table
col_w = 1.6
row_h = 0.75
x_start = 2.5
y_start = 6.5

# Column headers
for j, crit in enumerate(criteria):
    x = x_start + j * col_w + col_w / 2
    ax.text(x, y_start + 0.35, crit, ha="center", va="center", fontsize=8,
            fontweight="bold", color="white",
            bbox=dict(boxstyle="round,pad=0.2", facecolor=DARK_BLUE, edgecolor="none"))

# Row headers + cells
for i, method in enumerate(methods):
    y = y_start - (i + 1) * row_h
    bg = "#E8F5E9" if method == "pLIN" else (LIGHT_GRAY if i % 2 == 0 else "white")
    box = FancyBboxPatch((0.3, y - row_h/2 + 0.1), 14.0 - 0.6, row_h - 0.05,
                         boxstyle="round,pad=0.02",
                         facecolor=bg, edgecolor="#E0E0E0", linewidth=0.5)
    ax.add_patch(box)

    method_color = GREEN if method == "pLIN" else "#333333"
    method_weight = "bold" if method == "pLIN" else "normal"
    ax.text(1.4, y, method, ha="center", va="center", fontsize=10,
            fontweight=method_weight, color=method_color)

    for j, score in enumerate(scores[i]):
        x = x_start + j * col_w + col_w / 2
        if score == 1:
            symbol = "✓"
            color = GREEN
        elif score == 0.5:
            symbol = "~"
            color = ORANGE
        else:
            symbol = "✗"
            color = RED
        ax.text(x, y, symbol, ha="center", va="center", fontsize=14,
                fontweight="bold", color=color)

# Legend
ax.text(2, 0.3, "✓ = Full support", fontsize=9, color=GREEN, fontweight="bold")
ax.text(5, 0.3, "~ = Partial support", fontsize=9, color=ORANGE, fontweight="bold")
ax.text(8.5, 0.3, "✗ = Not supported", fontsize=9, color=RED, fontweight="bold")

for fmt in ["png", "pdf"]:
    fig10.savefig(os.path.join(FIG_DIR, f"Figure10_method_comparison.{fmt}"))
plt.close(fig10)
print("  Saved Figure10_method_comparison.png/pdf")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 11: v2.1 Upgrade Overview (Phase 1 + Phase 2)
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 11: v2.1 Upgrades ...")

fig11, ax = plt.subplots(1, 1, figsize=(14, 8))
ax.axis("off")
ax.set_xlim(0, 14)
ax.set_ylim(0, 8)

ax.text(7, 7.6, "Figure 11. pLIN v2.1 Feature Upgrades",
        ha="center", fontsize=14, fontweight="bold", color=DARK_BLUE)

# Phase 1
phase1_box = FancyBboxPatch((0.3, 4.5), 6.2, 2.8, boxstyle="round,pad=0.15",
                            facecolor="white", edgecolor=BLUE, linewidth=2.5)
ax.add_patch(phase1_box)
ax.text(3.4, 7.0, "Phase 1 — Immediate, High-Value", ha="center", fontsize=12,
        fontweight="bold", color=BLUE)

phase1 = [
    ("Adaptive Calibration", "Default ON — per-Inc\nthreshold calibration", BLUE),
    ("Length Warning", "Flags <5 kb plasmids\nwith noisy 4-mer profiles", ORANGE),
    ("Metadata Upload", "CSV/TSV with dates,\nlocation, patient data", GREEN),
    ("Mash ANI", "MinHash ANI estimation\nk=21, s=10,000", PURPLE),
]

for i, (title, desc, color) in enumerate(phase1):
    x = 0.6 + (i % 2) * 3.1
    y = 6.2 - (i // 2) * 1.3
    add_fancy_box(ax, x, y, 2.8, 1.0, color, f"{title}\n\n{desc}", fontsize=8)

# Phase 2
phase2_box = FancyBboxPatch((7.3, 4.5), 6.4, 2.8, boxstyle="round,pad=0.15",
                            facecolor="white", edgecolor=RED, linewidth=2.5)
ax.add_patch(phase2_box)
ax.text(10.5, 7.0, "Phase 2 — Advanced Genomic Resolution", ha="center", fontsize=12,
        fontweight="bold", color=RED)

phase2 = [
    ("FastANI", "True ANI computation\n--fragLen 1000 bp", BLUE),
    ("SNP Sub-typing", "minimap2 -cx asm5\nwithin L6 clusters", GREEN),
    ("Temporal Outbreak", "30-day window +\nAMR fingerprint match", RED),
]

for i, (title, desc, color) in enumerate(phase2):
    x = 7.6 + i * 2.1
    y = 5.6
    add_fancy_box(ax, x, y, 1.9, 1.2, color, f"{title}\n\n{desc}", fontsize=8)

# Impact summary at bottom
impact_box = FancyBboxPatch((0.3, 0.5), 13.4, 3.5, boxstyle="round,pad=0.15",
                            facecolor="#FAFAFA", edgecolor=GRAY, linewidth=1.5)
ax.add_patch(impact_box)
ax.text(7, 3.7, "Feature Impact Matrix", ha="center", fontsize=12,
        fontweight="bold", color=DARK_BLUE)

impacts = [
    ("Adaptive\nCalibration", "Built-in", "Better per-Inc\naccuracy", BLUE),
    ("Sequence\nWarning", "Built-in", "User\nconfidence", ORANGE),
    ("Metadata\nUpload", "Built-in", "Epi\ncontext", GREEN),
    ("Mash\nANI", "mash", "ANI\nvalidation", PURPLE),
    ("FastANI", "fastANI", "True\nANI", BLUE),
    ("SNP\nSub-type", "minimap2", "Outbreak\nresolution", GREEN),
    ("Temporal\nClusters", "Built-in", "Real-time\nsurveillance", RED),
]

for i, (feat, tool, impact, color) in enumerate(impacts):
    x = 0.6 + i * 1.9
    # Feature name
    add_fancy_box(ax, x, 2.5, 1.6, 0.8, color, feat, fontsize=7)
    # Tool
    tool_color = GREEN if tool == "Built-in" else GRAY
    ax.text(x + 0.8, 2.2, tool, ha="center", fontsize=7, fontweight="bold", color=tool_color)
    # Impact
    ax.text(x + 0.8, 1.0, impact, ha="center", fontsize=7, color="#333333")

ax.text(7, 1.7, "Required Tool:", ha="center", fontsize=8, color=GRAY, style="italic")

for fmt in ["png", "pdf"]:
    fig11.savefig(os.path.join(FIG_DIR, f"Figure11_v21_upgrades.{fmt}"))
plt.close(fig11)
print("  Saved Figure11_v21_upgrades.png/pdf")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 12: Inc Group Training Distribution + Classifier Performance
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 12: Inc Group Training & Performance ...")

fig12, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Training data distribution — read dynamically from pLIN assignments
import pandas as pd
PLIN_FILE = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")
try:
    plin_df = pd.read_csv(PLIN_FILE, sep="\t")
    inc_group_counts = plin_df["inc_type"].value_counts()
    inc_groups = [(name, count) for name, count in inc_group_counts.items()]
    total_plasmids = len(plin_df)
    n_unique_plin = plin_df["pLIN"].nunique()
    n_inc_groups = len(inc_groups)
except Exception:
    # Fallback to hardcoded values if file not found
    inc_groups = [
        ("IncFII", 4629), ("IncN", 1097), ("IncX1", 705), ("IncFIB", 97),
        ("ColRNAI", 91), ("IncF", 75), ("IncX3", 56), ("IncHI2", 36),
        ("IncI1", 27), ("IncI2", 25), ("IncX4", 24), ("IncR", 21),
        ("ColE", 19), ("IncC", 16), ("IncHI1", 16), ("IncFIC", 14),
        ("IncAC2", 14), ("IncA", 14), ("IncI", 11), ("IncFIBK", 11),
    ]
    total_plasmids = 6998
    n_unique_plin = 2454
    n_inc_groups = 20

names = [g[0] for g in inc_groups][::-1]
counts = [g[1] for g in inc_groups][::-1]
colors_bar = [BLUE if c > 100 else (ORANGE if c > 30 else RED) for c in counts]

ax1.barh(range(len(names)), counts, color=colors_bar, edgecolor="white", linewidth=0.5)
ax1.set_yticks(range(len(names)))
ax1.set_yticklabels(names, fontsize=8)
ax1.set_xlabel("Number of Training Sequences", fontsize=10)
ax1.set_title(f"A. Training Dataset Distribution ({n_inc_groups} Inc Groups)", fontsize=12,
              fontweight="bold", color=DARK_BLUE)

for i, (n, c) in enumerate(zip(names, counts)):
    ax1.text(c + max(counts) * 0.005, i, str(c), va="center", fontsize=7, color="#333333")

ax1.text(max(counts) * 0.65, 3, f"Total: {total_plasmids:,}\nsequences", fontsize=10,
         fontweight="bold", color=DARK_BLUE, ha="center",
         bbox=dict(boxstyle="round,pad=0.5", facecolor=LIGHT_BLUE, edgecolor=BLUE, alpha=0.8))

ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Panel B: Key performance metrics
ax2.axis("off")
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 10)

ax2.text(5, 9.3, "B. Key Performance Metrics", ha="center", fontsize=12,
         fontweight="bold", color=DARK_BLUE)

metrics = [
    ("Inc Detection Accuracy", "92.2%", "5-fold stratified CV, 20 groups", BLUE),
    ("Simpson's Diversity (D)", "0.979", "Strain-level discriminatory power", GREEN),
    ("XGBoost F1 Score", "0.903", "ML validation (nested CV)", ORANGE),
    ("Inc Concordance", "99.5%", "Composition vs known Inc groups", PURPLE),
    ("Unique pLIN Codes", f"{n_unique_plin:,}", "Strain-level (L6) resolution", RED),
    ("Training Plasmids", f"{total_plasmids:,}", f"Across {n_inc_groups} Inc groups", TEAL),
]

for i, (name, value, desc, color) in enumerate(metrics):
    y = 8.0 - i * 1.35
    # Value circle/box
    box = FancyBboxPatch((0.5, y - 0.4), 3.0, 0.9, boxstyle="round,pad=0.1",
                         facecolor=color, edgecolor="white", linewidth=1.5)
    ax2.add_patch(box)
    ax2.text(2.0, y + 0.1, value, ha="center", va="center", fontsize=16,
             fontweight="bold", color="white")

    # Name and description
    ax2.text(4.0, y + 0.15, name, va="center", fontsize=11, fontweight="bold", color=color)
    ax2.text(4.0, y - 0.2, desc, va="center", fontsize=9, color=GRAY)

fig12.tight_layout()
for fmt in ["png", "pdf"]:
    fig12.savefig(os.path.join(FIG_DIR, f"Figure12_inc_training_performance.{fmt}"))
plt.close(fig12)
print("  Saved Figure12_inc_training_performance.png/pdf")


# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("All manuscript figures generated successfully!")
print(f"Output directory: {FIG_DIR}")
print("=" * 60)
