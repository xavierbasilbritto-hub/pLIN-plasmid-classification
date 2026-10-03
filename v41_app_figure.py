#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Supplementary Figure 3: the pLIN app, from screenshots of a real run on the bundled Swiss
VIM-1 sample data (sample_data/swiss_vim1_outbreak; 8 FASTA files, 19 contigs; pLIN v4.1,
AMRFinderPlus and MOB-suite enabled). Screenshots: output/backbone_v41/figures/app_screens/.

Usage:
  python v41_app_figure.py
Output: output/backbone_v41/figures/SFig3.{pdf,png}
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(BASE_DIR, "output", "backbone_v41", "figures")
SCR = os.path.join(FIG, "app_screens")
MM = 1 / 25.4

# (file, x0, x1, y0, y1 in screenshot pixels, title)
PANELS = {
    "a": ("home.png", 360, 1380, 290, 1000, "Upload and options (v4.1 default)"),
    "b": ("Results.png", 372, 1360, 330, 880, "pLIN code per plasmid"),
    "c": ("Epidemiology.png", 372, 1360, 330, 2430, "Outbreak clusters and mobility"),
    "d": ("Cladogram.png", 372, 1360, 330, 1110, "Relatedness tree with codes"),
    "e": ("Export.png", 372, 1360, 330, 900, "Export with version stamp"),
}


def show(ax, key):
    f, x0, x1, y0, y1, title = PANELS[key]
    img = mpimg.imread(os.path.join(SCR, f))[y0:y1, x0:x1]
    ax.imshow(img, interpolation="lanczos")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#8a8985"); s.set_linewidth(0.6)
    ax.set_title(title, loc="left", fontsize=7, fontweight="bold", fontfamily="Arial")
    ax.text(-0.02, 1.0, key, transform=ax.transAxes, fontsize=9, fontweight="bold", fontfamily="Arial",
            ha="right", va="bottom")


def main():
    plt.rcParams.update({"pdf.fonttype": 42, "savefig.dpi": 300})
    fig = plt.figure(figsize=(180 * MM, 215 * MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.05], hspace=0.16, wspace=0.08)
    show(fig.add_subplot(gs[0, 0]), "a")
    show(fig.add_subplot(gs[0, 1]), "b")
    show(fig.add_subplot(gs[1:, 0]), "c")
    show(fig.add_subplot(gs[1, 1]), "d")
    show(fig.add_subplot(gs[2, 1]), "e")
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"SFig3.{ext}"), bbox_inches="tight", pad_inches=0.03)
    print("wrote SFig3")


if __name__ == "__main__":
    main()
