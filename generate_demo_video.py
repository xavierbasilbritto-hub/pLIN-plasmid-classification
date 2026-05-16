#!/usr/bin/env python3
"""
Generate a 5-7 minute pLIN demo video as an annotated slideshow MP4.
Uses matplotlib to render each frame, ffmpeg to stitch into video.

Output: output/pLIN_Demo_Video.mp4
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.gridspec import GridSpec
from PIL import Image
import subprocess
import glob
import shutil

# ── Paths ────────────────────────────────────────────────────────────────────
BASE    = os.path.dirname(os.path.abspath(__file__))
FRAMES  = os.path.join(BASE, "output", "video_frames")
OUT_MP4 = os.path.join(BASE, "output", "pLIN_Demo_Video.mp4")
LOGO    = os.path.join(BASE, "data", "plin_logo.png")
os.makedirs(FRAMES, exist_ok=True)

# ── Palette ──────────────────────────────────────────────────────────────────
C = {
    "dark_blue":  "#0D47A1",
    "med_blue":   "#1E88E5",
    "light_blue": "#BBDEFB",
    "accent":     "#E3F2FD",
    "navy":       "#0A2A6E",
    "eu_blue":    "#003F9A",
    "eu_yellow":  "#FFCC00",
    "white":      "#FFFFFF",
    "black":      "#000000",
    "dark_gray":  "#333333",
    "mid_gray":   "#757575",
    "light_gray": "#F5F5F5",
    "green":      "#2E7D32",
    "light_green":"#C8E6C9",
    "orange":     "#EF6C00",
    "red":        "#C62828",
    "teal":       "#00796B",
    "purple":     "#6A1B9A",
}

W, H = 1920, 1080
DPI  = 96
FIG_W = W / DPI
FIG_H = H / DPI

FRAME_IDX = [0]

def new_fig():
    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor=C["white"])
    return fig


def save_frame(fig, duration_s=4.0, fps=25):
    """Save fig as N duplicate frames for the given duration."""
    n = int(duration_s * fps)
    path = os.path.join(FRAMES, f"frame_{FRAME_IDX[0]:06d}.png")
    fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=0,
                facecolor=fig.get_facecolor())
    plt.close(fig)
    # Duplicate for duration
    for i in range(1, n):
        dest = os.path.join(FRAMES, f"frame_{FRAME_IDX[0]+i:06d}.png")
        shutil.copy2(path, dest)
    FRAME_IDX[0] += n
    print(f"  Frame {FRAME_IDX[0]-n:06d}  ({duration_s}s)")


def draw_header_bar(ax, title, subtitle=None, bar_color=None):
    """Draw a coloured header bar across the top of an axes."""
    bc = bar_color or C["dark_blue"]
    ax.add_patch(FancyBboxPatch((0, 0.88), 1, 0.12, boxstyle="square,pad=0",
                                fc=bc, ec="none", transform=ax.transAxes, zorder=3))
    ax.text(0.02, 0.94, title, transform=ax.transAxes,
            fontsize=28, fontweight="bold", color="white",
            va="center", ha="left", zorder=4)
    if subtitle:
        ax.text(0.02, 0.895, subtitle, transform=ax.transAxes,
                fontsize=14, color=C["light_blue"], style="italic",
                va="center", ha="left", zorder=4)
    # Yellow accent line
    ax.add_patch(FancyBboxPatch((0, 0.875), 1, 0.006, boxstyle="square,pad=0",
                                fc=C["eu_yellow"], ec="none",
                                transform=ax.transAxes, zorder=5))


def draw_footer(ax, left_text="pLIN v2.1 | UMCG | DRAIGON Consortium",
                right_text="Horizon Europe GA No. 101137383"):
    ax.add_patch(FancyBboxPatch((0, 0), 1, 0.045, boxstyle="square,pad=0",
                                fc=C["light_gray"], ec="none",
                                transform=ax.transAxes, zorder=3))
    ax.text(0.01, 0.022, left_text, transform=ax.transAxes,
            fontsize=10, color=C["mid_gray"], va="center", ha="left", zorder=4)
    ax.text(0.99, 0.022, right_text, transform=ax.transAxes,
            fontsize=10, color=C["mid_gray"], va="center", ha="right", zorder=4)


def pill(ax, x, y, w, h, text, fc, tc="white", fs=13, bold=False,
         transform=None, zorder=5, radius=0.015):
    tf = transform or ax.transAxes
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle=f"round,pad=0,rounding_size={radius}",
                                fc=fc, ec="none", transform=tf, zorder=zorder))
    ax.text(x + w/2, y + h/2, text, transform=tf,
            fontsize=fs, color=tc, va="center", ha="center",
            fontweight="bold" if bold else "normal", zorder=zorder+1)


def card(ax, x, y, w, h, title, value, detail,
         fc=None, tc="white", vc=None):
    fc = fc or C["dark_blue"]
    vc = vc or C["eu_yellow"]
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0,rounding_size=0.01",
                                fc=fc, ec="none", transform=ax.transAxes, zorder=3))
    ax.text(x + w/2, y + h - 0.025, title, transform=ax.transAxes,
            fontsize=12, color=tc, va="top", ha="center", zorder=4)
    ax.text(x + w/2, y + h/2 + 0.01, value, transform=ax.transAxes,
            fontsize=26, fontweight="bold", color=vc,
            va="center", ha="center", zorder=4)
    ax.text(x + w/2, y + 0.02, detail, transform=ax.transAxes,
            fontsize=10, color=tc, va="bottom", ha="center",
            style="italic", zorder=4)


# ════════════════════════════════════════════════════════════════════════════
# SCENE 1 — Title card (5 s)
# ════════════════════════════════════════════════════════════════════════════
print("Scene 1: Title")
fig = new_fig()
fig.patch.set_facecolor(C["navy"])
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.set_facecolor(C["navy"])

# EU funding bar
ax.add_patch(FancyBboxPatch((0, 0.95), 1, 0.05, boxstyle="square,pad=0",
                            fc=C["eu_blue"], ec="none", transform=ax.transAxes))
ax.text(0.5, 0.975, "★  Funded by the European Union  |  Horizon Europe  |  Grant Agreement No. 101137383  ★",
        transform=ax.transAxes, fontsize=13, color=C["eu_yellow"],
        va="center", ha="center")

# Logo-style title block
pill(ax, 0.25, 0.72, 0.5, 0.07, "DRAIGON ANNUAL CONSORTIUM MEETING",
     C["eu_yellow"], tc=C["navy"], fs=17, bold=True)

ax.text(0.5, 0.6, "pLIN", transform=ax.transAxes,
        fontsize=110, fontweight="bold", color="white",
        va="center", ha="center",
        path_effects=[pe.withStroke(linewidth=6, foreground=C["med_blue"])])

ax.text(0.5, 0.46, "Plasmid Lineage Identification Number",
        transform=ax.transAxes, fontsize=28, color=C["light_blue"],
        va="center", ha="center", style="italic")

ax.text(0.5, 0.38, "A Hierarchical Classification & AMR Surveillance Platform",
        transform=ax.transAxes, fontsize=18, color=C["eu_yellow"],
        va="center", ha="center")

# Stats preview
stats = [("79,305", "Plasmids"), ("28", "Inc/Rep\nGroups"),
         ("91.1%", "Classifier\nAccuracy"), ("0.985", "Simpson's D")]
for i, (n, l) in enumerate(stats):
    x0 = 0.1 + i * 0.2
    ax.add_patch(FancyBboxPatch((x0, 0.19), 0.17, 0.12,
                                boxstyle="round,pad=0,rounding_size=0.01",
                                fc=C["dark_blue"], ec=C["med_blue"],
                                linewidth=1.5, transform=ax.transAxes))
    ax.text(x0 + 0.085, 0.275, n, transform=ax.transAxes,
            fontsize=22, fontweight="bold", color=C["eu_yellow"],
            va="center", ha="center")
    ax.text(x0 + 0.085, 0.215, l, transform=ax.transAxes,
            fontsize=11, color="white", va="center", ha="center")

ax.text(0.5, 0.10, "Basil Britto Xavier, Anurag Kumar Bari, Bhanu Sinha, John W.A. Rossen",
        transform=ax.transAxes, fontsize=15, color="white",
        va="center", ha="center")
ax.text(0.5, 0.06, "Department of Medical Microbiology and Infection Prevention | UMCG, University of Groningen",
        transform=ax.transAxes, fontsize=12, color=C["light_blue"],
        va="center", ha="center")
save_frame(fig, duration_s=5)


# ════════════════════════════════════════════════════════════════════════════
# SCENE 2 — The Problem (5 s)
# ════════════════════════════════════════════════════════════════════════════
print("Scene 2: The Problem")
fig = new_fig()
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.set_facecolor(C["white"])
draw_header_bar(ax, "The Challenge: Plasmid Typing Gaps in AMR Surveillance",
                "Why current tools can't answer the outbreak question")
draw_footer(ax)

# Big AMR stat
ax.add_patch(FancyBboxPatch((0.03, 0.72), 0.94, 0.14,
                            boxstyle="round,pad=0,rounding_size=0.01",
                            fc=C["red"], ec="none", transform=ax.transAxes))
ax.text(0.5, 0.81, "4.95 million deaths associated with bacterial AMR in 2019",
        transform=ax.transAxes, fontsize=22, fontweight="bold",
        color="white", va="center", ha="center")
ax.text(0.5, 0.752, "Plasmid-mediated horizontal gene transfer is the #1 driver   |   Murray et al., Lancet 2022",
        transform=ax.transAxes, fontsize=13, color=C["light_gray"],
        va="center", ha="center", style="italic")

# Three problems
problems = [
    ("Limited Resolution", C["orange"],
     '"IncN" = >1,000 plasmids\nwith vastly different resistance\nprofiles — untraceable.'),
    ("Unstable Codes", C["red"],
     'MOB-suite / COPLA reassign\nidentifiers with each update —\nno longitudinal comparison.'),
    ("No AMR Integration", C["purple"],
     'Zero tools combine\nhierarchical typing AND\nresistance gene profiling.'),
]
for i, (title, col, desc) in enumerate(problems):
    x0 = 0.03 + i * 0.325
    ax.add_patch(FancyBboxPatch((x0, 0.33), 0.29, 0.36,
                                boxstyle="round,pad=0,rounding_size=0.01",
                                fc=C["light_gray"], ec=col, linewidth=2,
                                transform=ax.transAxes))
    ax.add_patch(FancyBboxPatch((x0, 0.685), 0.29, 0.007,
                                boxstyle="square,pad=0",
                                fc=col, ec="none", transform=ax.transAxes))
    ax.text(x0 + 0.145, 0.63, title, transform=ax.transAxes,
            fontsize=17, fontweight="bold", color=col,
            va="center", ha="center")
    ax.text(x0 + 0.145, 0.48, desc, transform=ax.transAxes,
            fontsize=13, color=C["dark_gray"],
            va="center", ha="center", linespacing=1.6)

# Solution callout
ax.add_patch(FancyBboxPatch((0.03, 0.1), 0.94, 0.2,
                            boxstyle="round,pad=0,rounding_size=0.01",
                            fc=C["accent"], ec=C["dark_blue"], linewidth=1.5,
                            transform=ax.transAxes))
ax.text(0.5, 0.25, "pLIN solves all three: permanent 6-level codes + integrated AMR profiling",
        transform=ax.transAxes, fontsize=20, fontweight="bold",
        color=C["dark_blue"], va="center", ha="center")
ax.text(0.5, 0.155, "Simpson's D = 0.985 vs 0.641 for Inc typing alone — 1.54× improvement in discriminatory power",
        transform=ax.transAxes, fontsize=14, color=C["dark_gray"],
        va="center", ha="center", style="italic")
save_frame(fig, duration_s=5)


# ════════════════════════════════════════════════════════════════════════════
# SCENE 3 — How pLIN Works (pipeline) (7 s)
# ════════════════════════════════════════════════════════════════════════════
print("Scene 3: Pipeline")
fig = new_fig()
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.set_facecolor(C["white"])
draw_header_bar(ax, "How pLIN Works: From Sequence to Hierarchical Code",
                "Tetranucleotide composition + single-linkage clustering at 6 ANI-calibrated thresholds")
draw_footer(ax)

# Pipeline steps
steps = [
    ("1", "Input\nFASTA",          C["dark_blue"],  "Any plasmid\nor contig"),
    ("2", "4-mer\nVectors",         C["med_blue"],   "256 features\nL2-normalised"),
    ("3", "Cosine\nDistances",      C["teal"],       "Scale-invariant\npairwise"),
    ("4", "Single-Linkage\nCluster",C["green"],      "6 thresholds\ncalibrated to ANI"),
    ("5", "Permanent\npLIN Code",   C["orange"],     "A.B.C.D.E.F\nidentifier"),
    ("6", "AMR\nReports",           C["red"],        "AMRFinderPlus\noutbreak risk"),
]
n = len(steps)
bw = 0.12
bh = 0.22
gap = (1.0 - n * bw) / (n + 1)
y0 = 0.56

for i, (num, title, col, desc) in enumerate(steps):
    x0 = gap + i * (bw + gap)
    ax.add_patch(FancyBboxPatch((x0, y0), bw, bh,
                                boxstyle="round,pad=0,rounding_size=0.01",
                                fc=col, ec="none", transform=ax.transAxes))
    ax.text(x0 + bw/2, y0 + bh - 0.03, num, transform=ax.transAxes,
            fontsize=26, fontweight="bold", color=C["eu_yellow"],
            va="top", ha="center")
    ax.text(x0 + bw/2, y0 + bh/2 - 0.01, title, transform=ax.transAxes,
            fontsize=12, fontweight="bold", color="white",
            va="center", ha="center", linespacing=1.5)
    ax.text(x0 + bw/2, y0 + 0.015, desc, transform=ax.transAxes,
            fontsize=9, color="white", va="bottom", ha="center",
            linespacing=1.4, style="italic")
    # Arrow
    if i < n - 1:
        ax.annotate("", xy=(x0 + bw + gap * 0.9, y0 + bh/2),
                    xytext=(x0 + bw + gap * 0.1, y0 + bh/2),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", color=C["mid_gray"],
                                   lw=2.5))

# 6-level table
levels = [
    ("L1", "≤ 0.150", "~85% ANI", "Superfamily",   C["dark_blue"]),
    ("L2", "≤ 0.100", "~90% ANI", "Major lineage",  C["med_blue"]),
    ("L3", "≤ 0.050", "~95% ANI", "Cluster",        C["teal"]),
    ("L4", "≤ 0.020", "~98% ANI", "Sublineage",     C["green"]),
    ("L5", "≤ 0.010", "~99% ANI", "Clone complex",  C["orange"]),
    ("L6", "≤ 0.001", "~99.9%",   "Strain/Outbreak",C["red"]),
]
# Table header
col_w = [0.06, 0.11, 0.11, 0.18]
col_x = [0.04, 0.10, 0.21, 0.32]
headers = ["Level", "Distance", "ANI", "Resolution"]
for j, (hdr, cx) in enumerate(zip(headers, col_x)):
    ax.add_patch(FancyBboxPatch((cx, 0.45), col_w[j], 0.065,
                                boxstyle="square,pad=0",
                                fc=C["dark_blue"], ec="none",
                                transform=ax.transAxes))
    ax.text(cx + col_w[j]/2, 0.483, hdr, transform=ax.transAxes,
            fontsize=12, fontweight="bold", color="white",
            va="center", ha="center")

for r, (lv, dist, ani, res, col) in enumerate(levels):
    bg = C["light_gray"] if r % 2 == 0 else C["white"]
    y_r = 0.45 - (r + 1) * 0.065
    for j, (val, cx) in enumerate(zip([lv, dist, ani, res], col_x)):
        ax.add_patch(FancyBboxPatch((cx, y_r), col_w[j], 0.065,
                                    boxstyle="square,pad=0",
                                    fc=bg, ec=C["light_gray"],
                                    linewidth=0.5, transform=ax.transAxes))
        bold = (j == 0)
        ax.text(cx + col_w[j]/2, y_r + 0.033, val, transform=ax.transAxes,
                fontsize=11, fontweight="bold" if bold else "normal",
                color=col if j == 0 else C["dark_gray"],
                va="center", ha="center")

# Example pLIN code box
ax.add_patch(FancyBboxPatch((0.55, 0.11), 0.41, 0.39,
                            boxstyle="round,pad=0,rounding_size=0.01",
                            fc=C["navy"], ec=C["med_blue"], linewidth=1.5,
                            transform=ax.transAxes))
ax.text(0.755, 0.46, "Example pLIN Code", transform=ax.transAxes,
        fontsize=14, fontweight="bold", color=C["eu_yellow"],
        va="center", ha="center")
ax.text(0.755, 0.38, "1.1.2.15.48.671", transform=ax.transAxes,
        fontsize=32, fontweight="bold", color=C["eu_yellow"],
        va="center", ha="center", family="monospace")
ax.text(0.755, 0.33, "IncN | blaKPC-2 carrier | n=90 | CRITICAL",
        transform=ax.transAxes, fontsize=11, color=C["light_blue"],
        va="center", ha="center", style="italic")
details = [
    ("L1.L2 = 1.1", "Superfamily / Major lineage"),
    ("L3 = 2",      "Species-level cluster"),
    ("L4 = 15",     "Sublineage"),
    ("L5 = 48",     "Clone complex"),
    ("L6 = 671",    "Strain / Outbreak level"),
]
for k, (code, meaning) in enumerate(details):
    yk = 0.27 - k * 0.040
    ax.text(0.58, yk, code, transform=ax.transAxes,
            fontsize=11, fontweight="bold", color=C["eu_yellow"],
            family="monospace", va="center")
    ax.text(0.72, yk, meaning, transform=ax.transAxes,
            fontsize=11, color=C["light_blue"], va="center")

# Key principle
ax.text(0.5, 0.075, "★  Codes are PERMANENT — never reassigned as the database grows  ★",
        transform=ax.transAxes, fontsize=15, fontweight="bold",
        color=C["dark_blue"], va="center", ha="center")
save_frame(fig, duration_s=7)


# ════════════════════════════════════════════════════════════════════════════
# SCENE 4 — GUI Demo mockup (8 s)
# ════════════════════════════════════════════════════════════════════════════
print("Scene 4: GUI walkthrough")
fig = new_fig()
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.set_facecolor(C["white"])
draw_header_bar(ax, "pLIN GUI: Interactive Streamlit Web Application",
                "Upload FASTA → Classify → AMR profile → Outbreak risk  |  http://localhost:8501")
draw_footer(ax)

# Browser chrome
ax.add_patch(FancyBboxPatch((0.02, 0.11), 0.96, 0.77,
                            boxstyle="round,pad=0,rounding_size=0.005",
                            fc=C["light_gray"], ec=C["mid_gray"], linewidth=1,
                            transform=ax.transAxes, zorder=2))
# Browser top bar
ax.add_patch(FancyBboxPatch((0.02, 0.85), 0.96, 0.04,
                            boxstyle="square,pad=0",
                            fc="#E0E0E0", ec="none", transform=ax.transAxes, zorder=3))
for cx, col in zip([0.04, 0.065, 0.09], ["#FF5F57","#FEBC2E","#28C840"]):
    ax.add_patch(plt.Circle((cx, 0.87), 0.008, fc=col, ec="none",
                            transform=ax.transAxes, zorder=4))
ax.add_patch(FancyBboxPatch((0.12, 0.856), 0.74, 0.028,
                            boxstyle="round,pad=0,rounding_size=0.003",
                            fc="white", ec="#CCCCCC", linewidth=0.5,
                            transform=ax.transAxes, zorder=4))
ax.text(0.49, 0.87, "localhost:8501  —  pLIN v2.1 | Plasmid Lineage Identification Number",
        transform=ax.transAxes, fontsize=10, color=C["dark_gray"],
        va="center", ha="center", zorder=5)

# Streamlit sidebar
ax.add_patch(FancyBboxPatch((0.02, 0.11), 0.22, 0.74,
                            boxstyle="square,pad=0",
                            fc="#F8F9FA", ec="#E0E0E0", linewidth=0.5,
                            transform=ax.transAxes, zorder=3))
ax.text(0.13, 0.81, "⚙  pLIN Settings", transform=ax.transAxes,
        fontsize=13, fontweight="bold", color=C["dark_blue"],
        va="center", ha="center", zorder=4)

sidebar_items = [
    ("Upload FASTA", C["dark_blue"]),
    ("Inc/Rep Classification", C["teal"]),
    ("AMR Detection", C["orange"]),
    ("Outbreak Detection", C["red"]),
    ("Assembly QC", C["green"]),
    ("Recombination", C["purple"]),
    ("Novel Groups", C["med_blue"]),
    ("MGE Boundaries", C["navy"]),
]
for k, (label, col) in enumerate(sidebar_items):
    yk = 0.76 - k * 0.07
    selected = (k == 0)
    bg = col if selected else "#EEEEEE"
    ax.add_patch(FancyBboxPatch((0.03, yk - 0.015), 0.20, 0.042,
                                boxstyle="round,pad=0,rounding_size=0.005",
                                fc=bg, ec="none", transform=ax.transAxes, zorder=4))
    ax.text(0.13, yk + 0.006, label, transform=ax.transAxes,
            fontsize=10, color="white" if selected else C["dark_gray"],
            fontweight="bold" if selected else "normal",
            va="center", ha="center", zorder=5)

# Main panel
# File upload widget
ax.add_patch(FancyBboxPatch((0.26, 0.72), 0.70, 0.10,
                            boxstyle="round,pad=0,rounding_size=0.005",
                            fc="white", ec=C["med_blue"], linewidth=1.5,
                            linestyle="dashed", transform=ax.transAxes, zorder=3))
ax.text(0.61, 0.795, "📁  Upload plasmid FASTA file",
        transform=ax.transAxes, fontsize=16, color=C["dark_blue"],
        va="center", ha="center", zorder=4)
ax.text(0.61, 0.755, "pKPC_demo.fasta  ✓  (47,293 bp | 1 sequence)",
        transform=ax.transAxes, fontsize=13, color=C["green"],
        va="center", ha="center", zorder=4)

# Classification output card
ax.add_patch(FancyBboxPatch((0.26, 0.52), 0.70, 0.175,
                            boxstyle="round,pad=0,rounding_size=0.005",
                            fc=C["accent"], ec=C["dark_blue"], linewidth=1,
                            transform=ax.transAxes, zorder=3))
ax.text(0.28, 0.675, "📊  Classification Result", transform=ax.transAxes,
        fontsize=14, fontweight="bold", color=C["dark_blue"],
        va="center", ha="left", zorder=4)

res_cols = [
    ("Inc/Rep Group", "IncN",            C["dark_blue"]),
    ("pLIN Code",     "1.1.2.15.48.671", C["teal"]),
    ("Confidence",    "99.2%",           C["green"]),
    ("Risk Tier",     "🔴 CRITICAL",      C["red"]),
]
for j, (lbl, val, col) in enumerate(res_cols):
    xj = 0.28 + j * 0.175
    ax.text(xj, 0.64, lbl, transform=ax.transAxes,
            fontsize=11, color=C["mid_gray"], va="center", ha="left")
    ax.text(xj, 0.575, val, transform=ax.transAxes,
            fontsize=16, fontweight="bold", color=col,
            va="center", ha="left")

# AMR section
ax.add_patch(FancyBboxPatch((0.26, 0.30), 0.33, 0.20,
                            boxstyle="round,pad=0,rounding_size=0.005",
                            fc="white", ec=C["orange"], linewidth=1,
                            transform=ax.transAxes, zorder=3))
ax.text(0.28, 0.475, "🧬  AMR Genes Detected", transform=ax.transAxes,
        fontsize=13, fontweight="bold", color=C["orange"],
        va="center", ha="left", zorder=4)
amr_genes = [
    ("blaKPC-2", "Carbapenemase", C["red"]),
    ("aac(6')-Ib", "Aminoglycoside", C["orange"]),
    ("qnrS1", "Quinolone (PMQR)", C["purple"]),
    ("blaTEM-1", "Beta-lactam", C["dark_blue"]),
]
for j, (gene, cls, col) in enumerate(amr_genes):
    yj = 0.43 - j * 0.034
    ax.text(0.285, yj, f"  •  {gene}", transform=ax.transAxes,
            fontsize=10, fontweight="bold", color=col, va="center")
    ax.text(0.42, yj, cls, transform=ax.transAxes,
            fontsize=10, color=C["mid_gray"], va="center", style="italic")

# Lineage info panel
ax.add_patch(FancyBboxPatch((0.61, 0.30), 0.35, 0.20,
                            boxstyle="round,pad=0,rounding_size=0.005",
                            fc=C["navy"], ec=C["med_blue"], linewidth=1,
                            transform=ax.transAxes, zorder=3))
ax.text(0.635, 0.475, "📍  Lineage Intelligence", transform=ax.transAxes,
        fontsize=13, fontweight="bold", color=C["eu_yellow"],
        va="center", ha="left", zorder=4)
lin_items = [
    "pLIN 671  ·  IncN  ·  n=90 members",
    "100% blaKPC-2 carriage",
    "Mean 13.2 AMR genes/plasmid",
    "Detected in 7 countries",
]
for j, item in enumerate(lin_items):
    yj = 0.44 - j * 0.033
    ax.text(0.635, yj, f"  ▸  {item}", transform=ax.transAxes,
            fontsize=10, color=C["light_blue"], va="center")

# Run button
pill(ax, 0.45, 0.175, 0.13, 0.055, "▶  Run Analysis",
     C["green"], tc="white", fs=14, bold=True)
ax.text(0.61, 0.20, "← Click to classify, profile AMR, and assign pLIN code",
        transform=ax.transAxes, fontsize=12, color=C["mid_gray"],
        va="center", ha="left", style="italic")
save_frame(fig, duration_s=8)


# ════════════════════════════════════════════════════════════════════════════
# SCENE 5 — AMR Landscape (6 s)
# ════════════════════════════════════════════════════════════════════════════
print("Scene 5: AMR Landscape")
fig = new_fig()
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.set_facecolor(C["white"])
draw_header_bar(ax, "AMR Landscape: 64,891 Resistance Gene Detections",
                "AMRFinderPlus v4.2.5 mapped to 79,305 plasmid lineages")
draw_footer(ax)

# 6 stat cards
stats6 = [
    ("64,891",  "Total gene\nhits",           C["dark_blue"]),
    ("83.1%",   "Plasmids\ncarrying AMR",     C["teal"]),
    ("1,635",   "Carbapenemase\ngenes",        C["red"]),
    ("1,804",   "ESBL\ngenes",                C["orange"]),
    ("204",     "Mobile colistin\n(mcr)",      C["purple"]),
    ("2,315",   "Quinolone\n(PMQR)",           C["green"]),
]
cw, ch = 0.14, 0.20
gap_c = (1.0 - 6 * cw) / 7
for i, (n, l, col) in enumerate(stats6):
    x0 = gap_c + i * (cw + gap_c)
    ax.add_patch(FancyBboxPatch((x0, 0.68), cw, ch,
                                boxstyle="round,pad=0,rounding_size=0.01",
                                fc=C["light_gray"], ec=col, linewidth=2,
                                transform=ax.transAxes))
    ax.add_patch(FancyBboxPatch((x0, 0.875), cw, 0.007,
                                boxstyle="square,pad=0",
                                fc=col, ec="none", transform=ax.transAxes))
    ax.text(x0 + cw/2, 0.80, n, transform=ax.transAxes,
            fontsize=22, fontweight="bold", color=col,
            va="center", ha="center")
    ax.text(x0 + cw/2, 0.715, l, transform=ax.transAxes,
            fontsize=10, color=C["dark_gray"],
            va="center", ha="center", linespacing=1.4)

# Carbapenemase bar chart
carba = [("blaKPC-2", 824), ("blaNDM-1", 228), ("blaKPC-3", 193),
         ("blaNDM-5", 89),  ("blaOXA-48", 61)]
max_v = 824
ax.text(0.04, 0.64, "Carbapenemase genes", transform=ax.transAxes,
        fontsize=14, fontweight="bold", color=C["dark_blue"])
for j, (gene, v) in enumerate(carba):
    yj = 0.595 - j * 0.065
    bw_j = 0.38 * v / max_v
    ax.add_patch(FancyBboxPatch((0.20, yj), bw_j, 0.045,
                                boxstyle="square,pad=0",
                                fc=C["red"], ec="none", transform=ax.transAxes))
    ax.text(0.185, yj + 0.022, gene, transform=ax.transAxes,
            fontsize=11, color=C["dark_gray"], va="center", ha="right",
            style="italic", fontweight="bold")
    ax.text(0.20 + bw_j + 0.01, yj + 0.022, str(v),
            transform=ax.transAxes, fontsize=11,
            color=C["red"], va="center", fontweight="bold")

# ESBL bar chart
esbl = [("blaCTX-M-15", 505), ("blaCTX-M-65", 319),
        ("blaSHV-12", 277),   ("mcr-1.1", 83)]
max_e = 505
ax.text(0.54, 0.64, "ESBL + Colistin (mcr)", transform=ax.transAxes,
        fontsize=14, fontweight="bold", color=C["dark_blue"])
for j, (gene, v) in enumerate(esbl):
    yj = 0.595 - j * 0.065
    bw_j = 0.34 * v / max_e
    col_b = C["orange"] if "mcr" not in gene else C["purple"]
    ax.add_patch(FancyBboxPatch((0.70, yj), bw_j, 0.045,
                                boxstyle="square,pad=0",
                                fc=col_b, ec="none", transform=ax.transAxes))
    ax.text(0.685, yj + 0.022, gene, transform=ax.transAxes,
            fontsize=11, color=C["dark_gray"], va="center", ha="right",
            style="italic", fontweight="bold")
    ax.text(0.70 + bw_j + 0.01, yj + 0.022, str(v),
            transform=ax.transAxes, fontsize=11,
            color=col_b, va="center", fontweight="bold")

# Bottom insight
ax.add_patch(FancyBboxPatch((0.03, 0.11), 0.94, 0.155,
                            boxstyle="round,pad=0,rounding_size=0.01",
                            fc=C["accent"], ec=C["dark_blue"], linewidth=1,
                            transform=ax.transAxes))
ax.text(0.5, 0.22, "Why this matters: every resistance gene is mapped to a permanent, hierarchical pLIN code",
        transform=ax.transAxes, fontsize=17, fontweight="bold",
        color=C["dark_blue"], va="center", ha="center")
ax.text(0.5, 0.155, "pLIN 671 (IncN, n=90): 100% blaKPC-2 carriage  |  pLIN 860 (multi-Inc, n=142): 44.4% mcr  |  mean 14.4 AMR genes/plasmid",
        transform=ax.transAxes, fontsize=13, color=C["dark_gray"],
        va="center", ha="center", style="italic")
save_frame(fig, duration_s=6)


# ════════════════════════════════════════════════════════════════════════════
# SCENE 6 — Outbreak Validation (6 s)
# ════════════════════════════════════════════════════════════════════════════
print("Scene 6: Outbreak Validation")
fig = new_fig()
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.set_facecolor(C["white"])
draw_header_bar(ax, "Outbreak Validation: 27 Studies, 13 Countries, 7 Mechanisms",
                "Independent cross-validation against published surveillance data")
draw_footer(ax)

# Stats strip
vstats = [
    ("74", "Outbreak\nplasmids", C["dark_blue"]),
    ("27", "Published\nstudies",  C["teal"]),
    ("13", "Countries\ntested",   C["orange"]),
    ("7",  "Resistance\nmechanisms", C["red"]),
    ("85.1%", "High-confidence\nclassification", C["green"]),
    ("9", "Intra-study\nclusters detected", C["purple"]),
]
vw = 0.14; vh = 0.18; vgap = (1.0 - 6*vw) / 7
for i, (n, l, col) in enumerate(vstats):
    x0 = vgap + i * (vw + vgap)
    ax.add_patch(FancyBboxPatch((x0, 0.70), vw, vh,
                                boxstyle="round,pad=0,rounding_size=0.01",
                                fc=col, ec="none", transform=ax.transAxes))
    ax.text(x0 + vw/2, 0.825, n, transform=ax.transAxes,
            fontsize=22, fontweight="bold", color=C["eu_yellow"],
            va="center", ha="center")
    ax.text(x0 + vw/2, 0.73, l, transform=ax.transAxes,
            fontsize=10, color="white",
            va="center", ha="center", linespacing=1.4)

# Two lineage boxes
lineages = [
    ("pLIN 671 — IncN",
     C["red"],
     [("Members",       "n = 90"),
      ("blaKPC-2 carriage", "100%"),
      ("Mean AMR genes","13.2"),
      ("Risk tier",     "🔴 CRITICAL")],
     "Dominant KPC-2 IncN lineage spanning\nmultiple countries — invisible under\nconventional 'IncN' typing alone."),
    ("pLIN 860 — Multi-Inc",
     C["orange"],
     [("Members",       "n = 142 (5 Inc groups)"),
      ("mcr carriage",  "44.4%"),
      ("Mean AMR genes","14.4"),
      ("Risk tier",     "🔴 CRITICAL")],
     "Cross-Inc lineage with mobile colistin\nresistance + carbapenemases —\ndemonstrates MDR-stacking on one backbone."),
]
for idx, (name, col, stats_l, note) in enumerate(lineages):
    x0 = 0.03 + idx * 0.485
    ax.add_patch(FancyBboxPatch((x0, 0.26), 0.45, 0.40,
                                boxstyle="round,pad=0,rounding_size=0.01",
                                fc=C["light_gray"], ec=col, linewidth=2,
                                transform=ax.transAxes))
    ax.add_patch(FancyBboxPatch((x0, 0.648), 0.45, 0.008,
                                boxstyle="square,pad=0",
                                fc=col, ec="none", transform=ax.transAxes))
    ax.text(x0 + 0.225, 0.62, name, transform=ax.transAxes,
            fontsize=18, fontweight="bold", color=col,
            va="center", ha="center")
    for j, (k, v) in enumerate(stats_l):
        yj = 0.575 - j * 0.065
        ax.text(x0 + 0.025, yj, f"{k}:", transform=ax.transAxes,
                fontsize=12, fontweight="bold", color=C["dark_gray"],
                va="center")
        ax.text(x0 + 0.225, yj, v, transform=ax.transAxes,
                fontsize=12, fontweight="bold", color=col,
                va="center", ha="center")
    ax.text(x0 + 0.225, 0.305, note, transform=ax.transAxes,
            fontsize=11, color=C["mid_gray"], va="center", ha="center",
            style="italic", linespacing=1.5)

# MLST + pLIN combined typing
ax.add_patch(FancyBboxPatch((0.03, 0.11), 0.94, 0.12,
                            boxstyle="round,pad=0,rounding_size=0.01",
                            fc=C["navy"], ec="none", transform=ax.transAxes))
ax.text(0.5, 0.185, "Combined MLST + pLIN typing distinguishes clonal spread vs. horizontal plasmid transfer",
        transform=ax.transAxes, fontsize=17, fontweight="bold",
        color=C["eu_yellow"], va="center", ha="center")
ax.text(0.5, 0.14, "Directly answers infection-control questions at DRAIGON clinical sites",
        transform=ax.transAxes, fontsize=13, color="white",
        va="center", ha="center", style="italic")
save_frame(fig, duration_s=6)


# ════════════════════════════════════════════════════════════════════════════
# SCENE 7 — 7 Analytical Modules (6 s)
# ════════════════════════════════════════════════════════════════════════════
print("Scene 7: Analytical Modules")
fig = new_fig()
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.set_facecolor(C["white"])
draw_header_bar(ax, "Seven Analytical Modules: Beyond Classification",
                "Quality control, novelty detection, recombination, MGE, evolutionary dynamics")
draw_footer(ax)

modules = [
    ("L3", "Assembly\nCompleteness",  "94.2%\ncorrect QC",    C["dark_blue"]),
    ("L4", "Database\nCoverage",      "Traffic-light\nnovelty", C["teal"]),
    ("L5", "Recombination\nDetect.",  "87.3%/94.1%\nsens/spec", C["orange"]),
    ("L6", "Novel Group\nDiscovery",  "18 putative\nnew groups", C["purple"]),
    ("L7", "Evolutionary\nRate",      "1.8–8.4×10⁻⁶\nsubs/yr",  C["green"]),
    ("L8", "Cluster\nStability",      "ARI >0.85\n20/28 groups",C["red"]),
    ("L10","MGE Boundary\nDetect.",   "94.7% IS\nsensitivity",  C["navy"]),
]

bw2, bh2 = 0.123, 0.48
gap2 = (1.0 - 7 * bw2) / 8

for i, (lid, name, perf, col) in enumerate(modules):
    x0 = gap2 + i * (bw2 + gap2)
    ax.add_patch(FancyBboxPatch((x0, 0.18), bw2, bh2,
                                boxstyle="round,pad=0,rounding_size=0.01",
                                fc=C["light_gray"], ec=col, linewidth=1.5,
                                transform=ax.transAxes))
    # Top colour block
    ax.add_patch(FancyBboxPatch((x0, 0.18 + bh2 - 0.085), bw2, 0.085,
                                boxstyle="square,pad=0",
                                fc=col, ec="none", transform=ax.transAxes))
    ax.text(x0 + bw2/2, 0.18 + bh2 - 0.042, lid, transform=ax.transAxes,
            fontsize=20, fontweight="bold", color=C["eu_yellow"],
            va="center", ha="center")
    ax.text(x0 + bw2/2, 0.18 + bh2/2 - 0.01, name, transform=ax.transAxes,
            fontsize=11, fontweight="bold", color=C["dark_gray"],
            va="center", ha="center", linespacing=1.5)
    # Performance badge
    ax.add_patch(FancyBboxPatch((x0 + 0.005, 0.19), bw2 - 0.01, 0.068,
                                boxstyle="round,pad=0,rounding_size=0.01",
                                fc=col, ec="none", transform=ax.transAxes))
    ax.text(x0 + bw2/2, 0.224, perf, transform=ax.transAxes,
            fontsize=10, fontweight="bold", color="white",
            va="center", ha="center", linespacing=1.4)

# DRAIGON relevance box
ax.add_patch(FancyBboxPatch((0.03, 0.085), 0.94, 0.075,
                            boxstyle="round,pad=0,rounding_size=0.01",
                            fc=C["accent"], ec=C["dark_blue"], linewidth=1,
                            transform=ax.transAxes))
ax.text(0.5, 0.127, "These modules give DRAIGON sites quality-controlled, clinically interpretable plasmid intelligence",
        transform=ax.transAxes, fontsize=15, fontweight="bold",
        color=C["dark_blue"], va="center", ha="center")
ax.text(0.5, 0.098, "from assembly QC to IS-element mapping — automatically flagging mosaic, novel, and recombinant plasmids",
        transform=ax.transAxes, fontsize=12, color=C["dark_gray"],
        va="center", ha="center", style="italic")
save_frame(fig, duration_s=6)


# ════════════════════════════════════════════════════════════════════════════
# SCENE 8 — 28 Inc/Rep Groups Coverage (5 s)
# ════════════════════════════════════════════════════════════════════════════
print("Scene 8: Database coverage")
fig = new_fig()
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.set_facecolor(C["white"])
draw_header_bar(ax, "28 Inc/Rep Groups: Broadest Coverage to Date",
                "Gram-negative · Gram-positive · WHO ESKAPE priority pathogens")
draw_footer(ax)

categories = [
    ("Gram-Negative\n(20 groups)", C["dark_blue"], [
        "IncFII (n=4,629)", "IncN (n=1,097)", "IncX1 (n=705)",
        "IncFIB, ColRNAI, IncF", "IncX3, IncHI2, IncI1, IncI2",
        "IncX4, IncR, ColE, IncC", "IncHI1, IncFIC, IncAC2",
        "IncA, IncI, IncFIBK",
    ]),
    ("Gram-Positive\n(4 groups)", C["teal"], [
        "repSA_large", "(S. aureus large plasmids)",
        "repSA_small", "(S. aureus small plasmids)",
        "repEF_conj", "(Enterococcus conjugative)",
        "repEF_res", "(Enterococcus resistance)",
    ]),
    ("Acinetobacter\nbaumannii (2)", C["orange"], [
        "repAci1",
        "(Small cryptic plasmids)",
        "repAci_large",
        "(Resistance / conjugative)",
    ]),
    ("Pseudomonas\naeruginosa (2)", C["red"], [
        "repPae_large",
        "(Megaplasmids / medium)",
        "repPae_small",
        "(Small plasmids)",
    ]),
]

cw3 = 0.22; ch3 = 0.54; gap3 = (1.0 - 4 * cw3) / 5
for i, (cat, col, items) in enumerate(categories):
    x0 = gap3 + i * (cw3 + gap3)
    ax.add_patch(FancyBboxPatch((x0, 0.2), cw3, ch3,
                                boxstyle="round,pad=0,rounding_size=0.01",
                                fc=C["light_gray"], ec=col, linewidth=2,
                                transform=ax.transAxes))
    ax.add_patch(FancyBboxPatch((x0, 0.2 + ch3 - 0.105), cw3, 0.105,
                                boxstyle="square,pad=0",
                                fc=col, ec="none", transform=ax.transAxes))
    ax.text(x0 + cw3/2, 0.2 + ch3 - 0.052, cat, transform=ax.transAxes,
            fontsize=13, fontweight="bold", color="white",
            va="center", ha="center", linespacing=1.5)
    for j, item in enumerate(items):
        yj = 0.2 + ch3 - 0.125 - j * 0.055
        ax.text(x0 + 0.015, yj, item, transform=ax.transAxes,
                fontsize=10.5, color=C["dark_gray"], va="center")

# Bottom stats
pill(ax, 0.08, 0.11, 0.18, 0.055, "79,305 plasmids", C["dark_blue"],
     fs=14, bold=True)
pill(ax, 0.30, 0.11, 0.18, 0.055, "57,886 unique codes", C["teal"],
     fs=14, bold=True)
pill(ax, 0.52, 0.11, 0.16, 0.055, "97.3% classified", C["green"],
     fs=14, bold=True)
pill(ax, 0.72, 0.11, 0.20, 0.055, "PLSDB 2025 + NCBI RefSeq", C["orange"],
     fs=13, bold=True)
save_frame(fig, duration_s=5)


# ════════════════════════════════════════════════════════════════════════════
# SCENE 9 — Deployment & Summary (5 s)
# ════════════════════════════════════════════════════════════════════════════
print("Scene 9: Deployment")
fig = new_fig()
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.set_facecolor(C["white"])
draw_header_bar(ax, "Deployment: Ready for DRAIGON Clinical Sites",
                "Open-source · Cross-platform · <30 minutes · No specialist bioinformatics required")
draw_footer(ax)

platforms = [
    ("macOS\n🍎", C["dark_blue"],  "./pLIN_macOS.sh"),
    ("Linux\n🐧",  C["teal"],      "./pLIN_Linux.sh"),
    ("Windows\n🪟",C["med_blue"], "pLIN_Windows.bat"),
]
for i, (name, col, cmd) in enumerate(platforms):
    x0 = 0.05 + i * 0.22
    ax.add_patch(FancyBboxPatch((x0, 0.62), 0.19, 0.24,
                                boxstyle="round,pad=0,rounding_size=0.01",
                                fc=col, ec="none", transform=ax.transAxes))
    ax.text(x0 + 0.095, 0.80, name, transform=ax.transAxes,
            fontsize=18, fontweight="bold", color="white",
            va="center", ha="center", linespacing=1.5)
    ax.add_patch(FancyBboxPatch((x0 + 0.01, 0.63), 0.17, 0.04,
                                boxstyle="round,pad=0,rounding_size=0.005",
                                fc=C["navy"], ec="none", transform=ax.transAxes))
    ax.text(x0 + 0.095, 0.65, cmd, transform=ax.transAxes,
            fontsize=11, color=C["eu_yellow"],
            va="center", ha="center", family="monospace")

# Workflow steps
ax.add_patch(FancyBboxPatch((0.72, 0.54), 0.26, 0.33,
                            boxstyle="round,pad=0,rounding_size=0.01",
                            fc=C["accent"], ec=C["dark_blue"], linewidth=1,
                            transform=ax.transAxes))
ax.text(0.85, 0.845, "CLI Workflow", transform=ax.transAxes,
        fontsize=14, fontweight="bold", color=C["dark_blue"],
        va="center", ha="center")
cli_steps = [
    "streamlit run plin_app.py",
    "→ Upload FASTA",
    "→ Auto-classify Inc/Rep",
    "→ Assign pLIN code",
    "→ Profile AMR",
    "→ Download report",
]
for j, step in enumerate(cli_steps):
    ax.text(0.74, 0.805 - j * 0.045, step, transform=ax.transAxes,
            fontsize=11, color=C["dark_blue"] if j == 0 else C["dark_gray"],
            fontweight="bold" if j == 0 else "normal",
            family="monospace" if j == 0 else "sans-serif")

# Summary checkboxes
checks = [
    ("✓", "Open-source (GPL-3.0) on GitHub",                      C["green"]),
    ("✓", "Runs on standard laptop / workstation (16 GB RAM)",      C["green"]),
    ("✓", "<30 minutes full pipeline: classify + AMR + risk tier",  C["green"]),
    ("✓", "Validated: 27 outbreak studies, 13 countries",           C["green"]),
    ("✓", "Permanent codes — longitudinal comparison across sites", C["green"]),
    ("✓", "79,305 plasmid reference database",                      C["green"]),
    ("→", "FastANI validation Gram-pos/ESKAPE groups (in progress)",C["orange"]),
    ("→", "Long-read / Nanopore integration (Year 3)",              C["orange"]),
]
ax.text(0.05, 0.53, "Deployment Checklist", transform=ax.transAxes,
        fontsize=16, fontweight="bold", color=C["dark_blue"])
for j, (icon, text, col) in enumerate(checks):
    yj = 0.49 - j * 0.052
    ax.text(0.05, yj, icon, transform=ax.transAxes,
            fontsize=14, fontweight="bold", color=col, va="center")
    ax.text(0.08, yj, text, transform=ax.transAxes,
            fontsize=12, color=C["dark_gray"], va="center")

# GitHub bar
ax.add_patch(FancyBboxPatch((0.03, 0.11), 0.94, 0.06,
                            boxstyle="round,pad=0,rounding_size=0.01",
                            fc=C["navy"], ec="none", transform=ax.transAxes))
ax.text(0.5, 0.14, "github.com/xavierbasilbritto-hub/pLIN-plasmid-classification  |  GPL-3.0  |  Free for DRAIGON partners",
        transform=ax.transAxes, fontsize=14, fontweight="bold",
        color=C["eu_yellow"], va="center", ha="center")
save_frame(fig, duration_s=5)


# ════════════════════════════════════════════════════════════════════════════
# SCENE 10 — Closing card (6 s)
# ════════════════════════════════════════════════════════════════════════════
print("Scene 10: Closing")
fig = new_fig()
fig.patch.set_facecolor(C["navy"])
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.set_facecolor(C["navy"])

ax.add_patch(FancyBboxPatch((0, 0.95), 1, 0.05, boxstyle="square,pad=0",
                            fc=C["eu_blue"], ec="none", transform=ax.transAxes))
ax.text(0.5, 0.975, "★  Funded by the European Union  |  Horizon Europe  |  Grant Agreement No. 101137383  ★",
        transform=ax.transAxes, fontsize=13, color=C["eu_yellow"],
        va="center", ha="center")

ax.text(0.5, 0.78, "Thank you", transform=ax.transAxes,
        fontsize=80, fontweight="bold", color="white",
        va="center", ha="center")
ax.text(0.5, 0.68, "Questions & Discussion",
        transform=ax.transAxes, fontsize=26, color=C["eu_yellow"],
        va="center", ha="center", style="italic")

# Summary pills
summary = [
    ("79,305 Plasmids",        C["dark_blue"]),
    ("28 Inc/Rep Groups",      C["teal"]),
    ("91.1% Accuracy",         C["green"]),
    ("64,891 AMR Hits",        C["orange"]),
    ("27 Outbreak Studies",    C["red"]),
    ("7 Analytical Modules",   C["purple"]),
]
tot_w = 6 * 0.13 + 5 * 0.02
start_x2 = (1 - tot_w) / 2
for i, (label, col) in enumerate(summary):
    x0 = start_x2 + i * (0.13 + 0.02)
    pill(ax, x0, 0.50, 0.13, 0.055, label, col, fs=12, bold=True)

ax.text(0.5, 0.39,
        "Basil Britto Xavier, Anurag Kumar Bari, Bhanu Sinha, John W.A. Rossen",
        transform=ax.transAxes, fontsize=17, fontweight="bold", color="white",
        va="center", ha="center")
ax.text(0.5, 0.34, "on behalf of the DRAIGON Consortium",
        transform=ax.transAxes, fontsize=14, color=C["light_blue"],
        va="center", ha="center", style="italic")
ax.text(0.5, 0.28, "Department of Medical Microbiology and Infection Prevention",
        transform=ax.transAxes, fontsize=13, color=C["light_blue"],
        va="center", ha="center")
ax.text(0.5, 0.24, "AGE Research Group | UMCG, University of Groningen, The Netherlands",
        transform=ax.transAxes, fontsize=13, color=C["light_blue"],
        va="center", ha="center")
ax.text(0.5, 0.175,
        "basilbritto.xavier@umcg.nl  |  github.com/xavierbasilbritto-hub/pLIN-plasmid-classification",
        transform=ax.transAxes, fontsize=14, fontweight="bold",
        color=C["eu_yellow"], va="center", ha="center")

ax.text(0.5, 0.095,
        "Views and opinions expressed are those of the authors only and do not necessarily reflect "
        "those of the European Union or HADEA.",
        transform=ax.transAxes, fontsize=10, color=C["mid_gray"],
        va="center", ha="center", style="italic")
save_frame(fig, duration_s=6)


# ════════════════════════════════════════════════════════════════════════════
# STITCH with ffmpeg
# ════════════════════════════════════════════════════════════════════════════
total_frames = FRAME_IDX[0]
print(f"\nTotal frames generated: {total_frames}")
print(f"Estimated duration: {total_frames / 25:.1f} seconds")

print("\nStitching video with ffmpeg...")
cmd = [
    "/opt/homebrew/bin/ffmpeg", "-y",
    "-framerate", "25",
    "-pattern_type", "glob",
    "-i", os.path.join(FRAMES, "frame_*.png"),
    "-c:v", "libx264",
    "-preset", "slow",
    "-crf", "18",
    "-pix_fmt", "yuv420p",
    "-vf", "scale=1920:1080",
    OUT_MP4,
]
result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode == 0:
    size_mb = os.path.getsize(OUT_MP4) / (1024 * 1024)
    print(f"\n✓  Video saved: {OUT_MP4}")
    print(f"   Size: {size_mb:.1f} MB  |  Duration: {total_frames/25:.1f}s  |  1920×1080 H.264")
else:
    print("ffmpeg stderr:", result.stderr[-800:])

# Clean up frames
shutil.rmtree(FRAMES)
print("Frames cleaned up.")
