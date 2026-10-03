#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Manuscript figures for pLIN v4.1. Every number is read from the result files under
output/backbone_v41 (the same files docs/PLIN_FACTS.json is built from); nothing is typed in.

Main figures
  Fig1  the code and how it is assigned
  Fig2  pre-registered confirmatory evaluation
  Fig3  stability, robustness, short reads and speed
  Fig4  hospital outbreaks (Swiss VIM-1, comparative study, leave-one-out, prospective simulation)
  Fig5  One Health clones
Supplementary figures
  SFig1 database build-up
  SFig2 core architecture (assignment rules)
  SFig4 all levels and all comparators on the confirmatory test set
  SFig5 comparative outbreak study and leave-one-out, all levels
  SFig6 prospective simulation, isolate by isolate
  SFig7 agreement with plasmid taxonomic units (PTUs)
(SFig3, the app, is made from screenshots by v41_app_figure.py.)

Usage:
  python v41_figures.py
Output: output/backbone_v41/figures/{Fig*,SFig*}.{pdf,png}
"""

import json
import os
import subprocess

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V41 = os.path.join(BASE_DIR, "output", "backbone_v41")
OUT = os.path.join(V41, "figures")

# Validated categorical slots (all-pairs, light surface): one colour per method family.
PLIN, MOB, PLING, MGE = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
LEGACY = "#8a8985"                                   # pLIN v4, the pre-registered baseline
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8985", "#e4e3df", "#ffffff"
RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]   # blue, light to dark
LEVEL_RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]       # ordinal L1..L6
MM = 1 / 25.4

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.dpi": 300, "figure.facecolor": SURF, "axes.facecolor": SURF,
})


# ----------------------------------------------------------------------------- helpers
def tsv(*p):
    return pd.read_csv(os.path.join(V41, *p), sep="\t")


def js(*p):
    return json.load(open(os.path.join(V41, *p)))


def style(ax, grid="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.5, zorder=0)
        ax.set_axisbelow(True)


def label(ax, letter, x=-0.12, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")


def box(ax, x, y, w, h, text, fc="#f4f3f0", ec=MUTED, fs=6.5, weight="normal", color=INK, ha="center"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.012",
                                fc=fc, ec=ec, lw=0.6, zorder=2))
    tx = x + w / 2 if ha == "center" else x + 0.012
    ax.text(tx, y + h / 2, text, ha=ha, va="center", fontsize=fs, fontweight=weight, color=color, zorder=3,
            linespacing=1.25)


def arrow(ax, x0, y0, x1, y1, text=None, side="above"):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.7, shrinkA=0, shrinkB=0, mutation_scale=7))
    if text:
        dy = 0.018 if side == "above" else -0.018
        ax.text((x0 + x1) / 2, (y0 + y1) / 2 + dy, text, ha="center", va="bottom" if dy > 0 else "top",
                fontsize=6, color=INK2)


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote", name)


def git_time(sha):
    return subprocess.run(["git", "show", "-s", "--format=%ad", "--date=format:%d %b %Y, %H:%M", sha],
                          cwd=BASE_DIR, capture_output=True, text=True).stdout.strip()


def metric(tm, method, truth, col="F1_w"):
    r = tm[(tm.method == method) & (tm.truth == truth)].iloc[0]
    return r[col] if col else r


LEVEL_TABLE = [  # (level, meaning, similarity, threshold, rule); fixed by the pre-registration
    ("L1", "backbone family", "shared protein families", "0.40", "founder"),
    ("L2", "backbone group", "shared protein families", "0.60", "founder"),
    ("L3", "shared backbone", "k-mer containment", "0.50", "founder"),
    ("L4", "backbone variant", "k-mer containment", "0.80", "founder"),
    ("L5", "lineage", "k-mer similarity", "0.80", "nearest neighbour"),
    ("L6", "outbreak clone", "k-mer similarity", "0.95", "nearest neighbour"),
]


# ----------------------------------------------------------------------------- Figure 1
def fig1():
    ex = tsv("onehealth", "cross_sector_examples.tsv")
    ex = ex[(ex.level == "L6")].sort_values("plasmids_with_gene", ascending=False).iloc[0]
    code = ex.lineage.split(".")
    fig = plt.figure(figsize=(180 * MM, 112 * MM))

    # a: anatomy of a code
    ax = fig.add_axes([0.0, 0.50, 1.0, 0.48]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    label(ax, "a", x=0.0, y=0.93)
    ax.text(0.03, 0.95, f"Example: the {ex.gene} outbreak clone, pLIN {ex.lineage}", fontsize=7, va="top", color=INK2)
    w, g, x0 = 0.148, 0.012, 0.03
    for k, (lv, meaning, sim, thr, rule) in enumerate(LEVEL_TABLE):
        x = x0 + k * (w + g)
        ax.add_patch(FancyBboxPatch((x, 0.60), w, 0.20, boxstyle="round,pad=0,rounding_size=0.02",
                                    fc=LEVEL_RAMP[k], ec="none"))
        ax.text(x + w / 2, 0.70, code[k], ha="center", va="center", fontsize=11, fontweight="bold",
                color="#ffffff" if k >= 1 else INK)
        ax.text(x + w / 2, 0.53, f"{lv}  {meaning}", ha="center", va="center", fontsize=7, fontweight="bold")
        ax.text(x + w / 2, 0.41, f"{sim}\n≥ {thr}", ha="center", va="center", fontsize=6.5, color=INK2,
                linespacing=1.3)
        ax.text(x + w / 2, 0.28, f"{rule} rule", ha="center", va="center", fontsize=6.5, color=INK2, style="italic")
        if k < 5:
            ax.text(x + w + g / 2, 0.70, ".", ha="center", va="center", fontsize=14, fontweight="bold", color=INK2)
    ax.annotate("", xy=(x0, 0.16), xytext=(x0 + 6 * w + 5 * g, 0.16),
                arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.7, mutation_scale=7))
    ax.text(x0, 0.08, "coarse: related backbones (protein content)", fontsize=6.5, color=INK2, ha="left", va="center")
    ax.text(x0 + 6 * w + 5 * g, 0.08, "fine: same lineage, outbreak clone (nucleotide k-mers)", fontsize=6.5,
            color=INK2, ha="right", va="center")

    # b: assignment workflow
    ax = fig.add_axes([0.0, 0.0, 1.0, 0.50]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    label(ax, "b", x=0.0, y=0.90)
    h = 0.22
    box(ax, 0.03, 0.50, 0.14, h, "Plasmid sequence\n(complete or\nassembled contigs)")
    box(ax, 0.22, 0.50, 0.15, h, "Sequence already\nin the database?\n(exact hash)")
    box(ax, 0.22, 0.08, 0.15, 0.20, "Return its published\ncode (instant)", fc="#e8f1fc", ec=PLIN)
    box(ax, 0.42, 0.62, 0.18, 0.20, "Genes (pyrodigal)\nprotein families: exact,\nthen MMseqs2 search")
    box(ax, 0.42, 0.36, 0.18, 0.20, "k-mer sketch\nFracMinHash, k = 21")
    box(ax, 0.65, 0.50, 0.15, h, "Nearest database\nplasmid: k-mer\nsimilarity \u2265 0.80?")
    box(ax, 0.85, 0.62, 0.13, 0.20, "Copy its L1 to L5;\nL6 if \u2265 0.95", fc="#e8f1fc", ec=PLIN)
    box(ax, 0.85, 0.28, 0.13, 0.24, "Founder rule for\nL1 to L4;\nnew L5 and L6", fc="#e8f1fc", ec=PLIN)
    arrow(ax, 0.17, 0.61, 0.22, 0.61)
    arrow(ax, 0.295, 0.50, 0.295, 0.28)
    ax.text(0.30, 0.39, "yes", fontsize=6, color=INK2, ha="left")
    arrow(ax, 0.37, 0.64, 0.42, 0.70); arrow(ax, 0.37, 0.58, 0.42, 0.47)
    ax.text(0.385, 0.61, "no", fontsize=6, color=INK2, ha="center", va="center")
    arrow(ax, 0.60, 0.72, 0.65, 0.64); arrow(ax, 0.60, 0.46, 0.65, 0.56)
    arrow(ax, 0.80, 0.66, 0.85, 0.72); arrow(ax, 0.80, 0.56, 0.85, 0.42)
    ax.text(0.825, 0.72, "yes", fontsize=6, color=INK2, ha="center"); ax.text(0.825, 0.45, "no", fontsize=6, color=INK2, ha="center")
    ax.text(0.42, 0.18, "Database codes never change: new plasmids are only added,\nso a published name stays valid in every later release.",
            fontsize=6.5, color=INK2, va="center")
    save(fig, "Fig1")


# ----------------------------------------------------------------------------- Figure 2
def forest(ax, rows, title, ref_name, diff):
    ys = np.arange(len(rows))[::-1]
    for y, (name, col, f, lo, hi) in zip(ys, rows):
        ax.plot([lo, hi], [y, y], color=col, lw=1.4, solid_capstyle="round", zorder=3)
        ax.plot(f, y, "o", ms=5.5, color=col, mec=SURF, mew=1.0, zorder=4)
        ax.text(hi + 0.02, y, f"{f:.2f}", va="center", fontsize=6.5, color=INK)
    ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlim(0, 1.05); ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel("Weighted pairwise F1 (95% CI)")
    ax.set_title(title, loc="left", fontweight="bold")
    style(ax, grid="x")
    ax.text(0.0, -0.25, f"Difference vs {ref_name}: {diff['difference']:+.2f} (95% CI {diff['lo']:.2f} to {diff['hi']:.2f})",
            transform=ax.transAxes, fontsize=6.5, color=INK2)


def fig2():
    tm, ep = tsv("confirm", "test_metrics.tsv"), js("confirm", "endpoints.json")
    best = lambda prefix, truth: max((m for m in tm.method.unique() if m.startswith(prefix)),
                                     key=lambda m: metric(tm, m, truth))
    fig = plt.figure(figsize=(180 * MM, 118 * MM))

    # a: study design
    ax = fig.add_axes([0.0, 0.76, 1.0, 0.22]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    label(ax, "a", x=0.0, y=0.92)
    t_pre = git_time("6e29113")
    n_test = len(tsv("confirm", "test_plasmids.tsv"))
    box(ax, 0.03, 0.10, 0.20, 0.70, "Development data only\n(benchmark, v4 evaluation,\noutbreak and VIM-1 sets);\nthresholds chosen here")
    box(ax, 0.28, 0.10, 0.20, 0.70, f"Design and analysis\nfrozen, pre-registered\n{t_pre}\n(commit 6e29113)", fc="#e8f1fc", ec=PLIN)
    box(ax, 0.53, 0.10, 0.20, 0.70, f"{n_test:,} fresh plasmids\nnever used in design;\n{ep['test_pairs']:,} pairs aligned (BLAST)")
    box(ax, 0.78, 0.10, 0.20, 0.70, f"Truth from alignment:\n{ep['same_lineage_pairs']} same-lineage pairs,\n{ep['related_backbone_pairs']:,} related-backbone pairs;\nevaluated once")
    for x in (0.23, 0.48, 0.73):
        arrow(ax, x, 0.45, x + 0.05, 0.45)

    # b, c: primary endpoints
    lin = [("pLIN v4.1 L5", PLIN, "pLIN v4.1 L5"), ("MOB-suite secondary", MOB, "MOB-suite secondary"),
           (f"pling ({best('pling', 'same_lineage').split()[-1]})", PLING, best("pling", "same_lineage")),
           ("mge-cluster", MGE, "mge-cluster"), (f"pLIN v4 {best('pLIN v4 ', 'same_lineage').split()[-1]}", LEGACY,
                                                  best("pLIN v4 ", "same_lineage"))]
    bb = [("pLIN v4.1 L1", PLIN, "pLIN v4.1 L1"), ("MOB-suite primary", MOB, "MOB-suite primary"),
          (f"pling ({best('pling', 'related_backbone').split()[-1]})", PLING, best("pling", "related_backbone")),
          ("mge-cluster", MGE, "mge-cluster"), ("pLIN v4 L1", LEGACY, "pLIN v4 L1")]
    rows = lambda spec, truth: [(n, c, *metric(tm, m, truth, None)[["F1_w", "F1_w_lo", "F1_w_hi"]]) for n, c, m in spec]
    ax = fig.add_axes([0.13, 0.15, 0.33, 0.48]); label(ax, "b", x=-0.36)
    forest(ax, rows(lin, "same_lineage"), "Same lineage (endpoint 2)", "MOB-suite secondary",
           ep["E2_lineage_v41L5_vs_MOBsecondary"])
    ax = fig.add_axes([0.63, 0.15, 0.33, 0.48]); label(ax, "c", x=-0.36)
    forest(ax, rows(bb, "related_backbone"), "Related backbone (endpoint 1)", "MOB-suite primary",
           ep["E1_backbone_v41L1_vs_MOBprimary"])
    save(fig, "Fig2")


# ----------------------------------------------------------------------------- Figure 3
def grouped(ax, groups, series, colors, names, ylabel, ns=None):
    x = np.arange(len(groups)); w = 0.36
    for k, (vals, col, nm) in enumerate(zip(series, colors, names)):
        xs = x + (k - 0.5) * w
        ax.bar(xs, vals, w - 0.04, color=col, label=nm, zorder=3)
        for xi, v in zip(xs, vals):
            ax.text(xi, v + 1.5, f"{v:g}", ha="center", va="bottom", fontsize=5.5, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g}\n(n = {n})" for g, n in zip(groups, ns)] if ns is not None else groups)
    ax.set_ylim(0, 115); ax.set_yticks([0, 25, 50, 75, 100]); ax.set_ylabel(ylabel)
    style(ax)


def fig3():
    rc, sp = js("confirm", "release_checks.json"), js("speed_v41.json")
    fr, vm = tsv("shortread", "fragmentation_metrics.tsv"), tsv("shortread", "vim_metrics.tsv")
    fig = plt.figure(figsize=(180 * MM, 120 * MM))

    # a: stability and re-query (headline numbers)
    ax = fig.add_axes([0.0, 0.56, 0.27, 0.40]); ax.axis("off"); label(ax, "a", x=0.0, y=1.0)
    runs = rc["R1_stability"]["runs"]
    tiles = [(f"{100 * rc['R1_stability']['min_retained_all_levels']:.0f}%",
              f"of database codes unchanged at every\nlevel as the database grew\n({len(runs)} runs: 25, 50, 75% snapshots x 5 seeds)"),
             (f"{int(rc['R2_reproduction']['identical'] * rc['R2_reproduction']['n'])}/{rc['R2_reproduction']['n']}",
              "fresh database plasmids re-typed from\nsequence returned their own code")]
    for k, (big, small) in enumerate(tiles):
        y = 0.80 - k * 0.48
        ax.text(0.04, y, big, fontsize=20, fontweight="bold", color=INK, va="center")
        ax.text(0.04, y - 0.18, small, fontsize=6.5, color=INK2, va="center", linespacing=1.3)

    # b: point mutations
    ax = fig.add_axes([0.36, 0.60, 0.26, 0.34]); label(ax, "b", x=-0.25)
    rob = rc["R3_robustness"]; lv = np.arange(1, 7)
    lab_y = {}
    for (rate, col, nm) in [("0.0001", RAMP[2], "0.01%"), ("0.001", RAMP[4], "0.1%"), ("0.01", RAMP[6], "1%")]:
        v = 100 * np.array(rob[rate])
        ax.plot(lv, v, "-o", color=col, lw=1.5, ms=4, mec=SURF, mew=0.8, zorder=3)
        lab_y[nm] = v[-1]
    ax.text(6.15, lab_y["0.1%"] + 4, f"0.1%: {lab_y['0.1%']:.0f}", va="center", fontsize=6, color=INK)
    ax.text(6.15, lab_y["0.01%"] - 5, f"0.01%: {lab_y['0.01%']:.0f}", va="center", fontsize=6, color=INK)
    ax.text(6.15, lab_y["1%"], "1%", va="center", fontsize=6, color=INK)
    ax.axhline(95, color=MUTED, lw=0.7, ls=(0, (3, 2)), zorder=2)
    ax.text(1, 50, "dashed: release criterion\nfor L1 to L5 (95%)", fontsize=5.8, color=INK2, va="top")
    ax.set_xticks(lv); ax.set_xticklabels([f"L{i}" for i in lv]); ax.set_xlim(0.8, 7.4)
    ax.set_ylim(-3, 105); ax.set_ylabel("Plasmids keeping their code (%)")
    ax.set_title("Random substitutions (100 plasmids)", loc="left", fontweight="bold")
    style(ax)

    # c: speed
    ax = fig.add_axes([0.74, 0.60, 0.24, 0.34]); label(ax, "c", x=-0.42)
    names = ["pLIN, plasmid\nalready in database", "pLIN, related\nto a database plasmid", "pLIN, divergent\n(new backbone)",
             "MOB-suite\n(mob_typer)"]
    vals = [sp["known"]["seconds_per_plasmid"], sp["related"]["seconds_per_plasmid"],
            sp["divergent"]["seconds_per_plasmid"], sp["MOB-suite_seconds_per_plasmid_8_threads"]]
    ys = np.arange(4)[::-1]
    ax.barh(ys, vals, 0.6, color=[PLIN, PLIN, PLIN, MOB], zorder=3)
    for y, v in zip(ys, vals):
        ax.text(v + 0.03, y, "< 0.01 s" if v < 0.01 else f"{v:.2f} s", va="center", fontsize=6.5)
    ax.set_yticks(ys); ax.set_yticklabels(names, fontsize=6)
    ax.set_xlim(0, 1.6); ax.set_xlabel(f"Seconds per plasmid ({sp['threads']} threads)")
    ax.set_title("Typing speed (100 plasmids each)", loc="left", fontweight="bold")
    style(ax, grid="x")

    # d: simulated Illumina reads, assembled with SPAdes
    ax = fig.add_axes([0.07, 0.09, 0.45, 0.33]); label(ax, "d", x=-0.11)
    f = fr[fr.group != "all"].set_index("group").loc[["1 contig", "2-5", "6-20", ">20", "not assembled"]]
    grouped(ax, ["1 contig", "2 to 5", "6 to 20", "> 20", "not\nassembled"], [f.L5_kept_pct.values, f.L6_kept_pct.values],
            [PLIN, RAMP[6]], ["L1 to L5", "L6"], "Same code as the complete\nsequence (%)", ns=f.plasmids.values)
    ax.set_xlabel("Contigs per plasmid after short-read assembly")
    ax.set_title(f"Simulated short reads, {int(fr[fr.group == 'all'].plasmids.iloc[0])} plasmids (SPAdes)",
                 loc="left", fontweight="bold")
    assert (f[["L1_kept_pct", "L2_kept_pct", "L3_kept_pct", "L4_kept_pct"]].values.T == f.L5_kept_pct.values).all()
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.80, 0.80))

    # e: Swiss VIM-1 genomes, simulated and real Illumina reads, plasmid bins from MOB-recon
    ax = fig.add_axes([0.62, 0.09, 0.36, 0.33]); label(ax, "e", x=-0.16)
    vr = tsv("shortread", "vim_real_metrics.tsv").set_index("set")
    v = vm.set_index("set")
    good, bad = v.index[v.index.str.startswith(">=")][0], v.index[v.index.str.startswith("<")][0]
    groups = [(v, good, "one bin,\nsimulated"), (vr, good, "one bin,\nreal reads"), (v, bad, "split or merged,\nsimulated"),
              (vr, bad, "split or merged,\nreal reads")]
    grouped(ax, [g for _, _, g in groups], [[t.loc[k, "L5_kept_pct"] for t, k, _ in groups], [t.loc[k, "L6_kept_pct"] for t, k, _ in groups]],
            [PLIN, RAMP[6]], ["L1 to L5", "L6"], "Same code as the long-read\nplasmid (%)", ns=[int(t.loc[k, "plasmids"]) for t, k, _ in groups])
    ax.tick_params(axis="x", labelsize=5.6)
    ax.set_xlabel("MOB-recon binning (\u2265 90% of k-mers in one bin) and read source")
    ax.set_title("Swiss VIM-1 genomes: Illumina reads, SPAdes, MOB-recon", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="upper right")
    save(fig, "Fig3")


# ----------------------------------------------------------------------------- Figure 4
def tradeoff(ax, df, xcol, ycol, xlab, ylab, log=False, plin_levels=(1, 3, 5, 6), offsets=None):
    offsets = offsets or {}
    groups = [("pLIN v4.1", PLIN), ("MOB-suite", MOB), ("pling", PLING), ("mge-cluster", MGE)]
    for prefix, col in groups:
        d = df[df.method.str.startswith(prefix)]
        if d.empty:
            continue
        if prefix == "pLIN v4.1":
            d = d[d.method.isin([f"pLIN v4.1 L{k}" for k in plin_levels])]
            ax.plot(d[xcol], d[ycol], "-", color=col, lw=1.0, zorder=2)
        ax.scatter(d[xcol], d[ycol], s=34, color=col, edgecolor=SURF, lw=1.0, zorder=3)
        for _, r in d.iterrows():
            nm = r.method.replace("pLIN v4.1 ", "pLIN ").replace("MOB-suite ", "MOB ")
            dx, dy, ha = offsets.get(r.method, (4, 0, "left"))
            ax.annotate(nm, (r[xcol], r[ycol]), xytext=(dx, dy), textcoords="offset points", fontsize=6,
                        va="center", ha=ha, color=INK)
    if log:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel(xlab); ax.set_ylabel(ylab)
    style(ax, grid="both")


def fig4():
    vp, vc = tsv("case_studies", "swiss_vim1_pairs.tsv"), tsv("case_studies", "swiss_vim1_codes.tsv")
    oc, loo, pr = tsv("case_studies", "outbreak_comparators.tsv"), tsv("case_studies", "loo_metrics.tsv"), \
        tsv("case_studies", "prospective_metrics.tsv")
    fig = plt.figure(figsize=(180 * MM, 150 * MM))

    # a: Swiss VIM-1, levels shared vs alignment
    ax = fig.add_axes([0.04, 0.53, 0.36, 0.42]); label(ax, "a", x=-0.08, y=1.03)
    ids = sorted(set(vp.contig_shorter) | set(vp.contig_longer), key=lambda s: int(s.split("VIM")[1].split("_")[0]))
    short = {i: "VIM" + i.split("VIM")[1].split("_")[0] for i in ids}
    n = len(ids); pos = {i: k for k, i in enumerate(ids)}
    agree = 0
    for _, r in vp.iterrows():
        a, b = sorted([pos[r.contig_shorter], pos[r.contig_longer]])
        lv = int(r.levels_shared_v41)
        ax.add_patch(Rectangle((a, n - 1 - b), 1, 1, fc=RAMP[lv] if lv else "#f4f3f0", ec=SURF, lw=1.5, zorder=2))
        ax.text(a + 0.5, n - 1 - b + 0.5, str(lv), ha="center", va="center", fontsize=6.5,
                color="#ffffff" if lv >= 4 else INK, zorder=4)
        ok = (lv >= 5) == bool(r.same_lineage_by_alignment)
        agree += ok
        if r.same_lineage_by_alignment:
            ax.add_patch(Rectangle((a + 0.06, n - 1 - b + 0.06), 0.88, 0.88, fc="none", ec=INK, lw=1.1, zorder=3))
        if not ok:
            ax.text(a + 0.86, n - 1 - b + 0.84, "*", ha="center", va="center", fontsize=8, color=INK, zorder=5)
    ax.set_xlim(0, n - 1); ax.set_ylim(0, n - 1); ax.set_aspect("equal")
    ax.set_xticks(np.arange(n - 1) + 0.5); ax.set_xticklabels([short[i] for i in ids[:-1]], rotation=90)
    ax.set_yticks(np.arange(n - 1) + 0.5); ax.set_yticklabels([short[i] for i in ids[1:]][::-1])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title(f"Swiss blaVIM-1 outbreak plasmids:\n{agree} of {len(vp)} pairs agree with alignment",
                 loc="left", fontweight="bold", fontsize=6.5)
    ax.text(n - 1, n - 1.6, "Number = pLIN levels shared\nOutlined = same lineage by alignment\n(aligned fraction ≥ 0.80,"
            " identity ≥ 99%)\n* = disagreement", ha="right", va="top", fontsize=6, color=INK2, linespacing=1.3)

    # b: comparative outbreak study
    ax = fig.add_axes([0.57, 0.57, 0.40, 0.37]); label(ax, "b", x=-0.19)
    d = oc[oc.method != "pling community"].copy()
    d["different_study_together_pct"] = d.different_study_together_pct.clip(lower=0.1)
    tradeoff(ax, d, "same_study_together_pct", "different_study_together_pct",
             "Same-outbreak pairs grouped together (%)", "Pairs from different outbreaks\ngrouped together (%, log)",
             log=True, offsets={"pLIN v4.1 L6": (-4, 0, "right"), "pLIN v4.1 L5": (4, -5, "left"),
                                "MOB-suite secondary": (-5, -1, "right"), "pLIN v4.1 L3": (4, -2, "left"),
                                "mge-cluster": (-4, 0, "right")})
    ax.set_xlim(30, 105)
    ax.set_title(f"Published outbreaks: {int(oc.same_study_pairs.iloc[0] + oc.different_study_pairs.iloc[0]):,} pairs, "
                 f"74 plasmids, 27 studies", loc="left", fontweight="bold")

    # c: leave-one-out
    ax = fig.add_axes([0.07, 0.08, 0.36, 0.33]); label(ax, "c", x=-0.17)
    lv = np.arange(1, 7)
    for col, c, nm in [("recovery_pct", PLIN, "own code recovered"), ("same_study_pairs_pct", RAMP[6], "same-outbreak pairs together")]:
        ax.plot(lv, loo[col], "-o", color=c, lw=1.5, ms=4, mec=SURF, mew=0.8, zorder=3)
    ax.text(6.15, loo.recovery_pct.iloc[-1], "own code\nrecovered", fontsize=6, va="center")
    ax.text(6.15, loo.same_study_pairs_pct.iloc[-1] - 6, "same-outbreak\npairs together", fontsize=6, va="center")
    for k in (4,):
        ax.annotate(f"{loo.recovery_pct.iloc[k]:.1f}%", (lv[k], loo.recovery_pct.iloc[k]), xytext=(0, 6),
                    textcoords="offset points", ha="center", fontsize=6)
    ax.set_xticks(lv); ax.set_xticklabels([f"L{i}" for i in lv]); ax.set_xlim(0.7, 7.3)
    ax.set_ylim(0, 105); ax.set_ylabel("%")
    ax.set_title(f"Leave-one-out: {int(loo.plasmids_in_database.iloc[0])} outbreak plasmids,\n"
                 "each removed from the database and re-typed", loc="left", fontweight="bold")
    style(ax)

    # d: prospective simulation
    ax = fig.add_axes([0.57, 0.08, 0.40, 0.33]); label(ax, "d", x=-0.19)
    tradeoff(ax, pr, "early_warning_pct", "false_alarm_pct", "Early warning: isolate linked to an earlier\n"
             "isolate of its own outbreak (%)", "False alarm: linked to an earlier isolate\nof another outbreak (%)",
             offsets={"pLIN v4.1 L6": (-4, 0, "right"), "pLIN v4.1 L5": (4, -5, "left"), "pLIN v4.1 L3": (4, -7, "left"),
                      "MOB-suite primary": (-4, 0, "right"), "MOB-suite secondary": (-4, 6, "right")})
    ax.text(0.02, 0.97, "MOB-suite clusters come from its fixed reference\ndatabase, which includes plasmids deposited later",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.8, color=INK2)
    ax.set_xlim(35, 80); ax.set_ylim(0, 50)
    ax.set_title(f"Prospective simulation: {int(pr.isolates_with_earlier_outbreak_isolate.iloc[0])} isolates typed "
                 "in deposition order", loc="left", fontweight="bold")
    save(fig, "Fig4")


# ----------------------------------------------------------------------------- Figure 5
SECTORS = ["human", "animal", "food", "environment", "plant"]


def fig5():
    oh, ex = js("onehealth", "summary.json"), tsv("onehealth", "cross_sector_examples.tsv")
    fig = plt.figure(figsize=(180 * MM, 95 * MM))

    # a: lineages crossing sectors
    ax = fig.add_axes([0.06, 0.16, 0.25, 0.66]); label(ax, "a", x=-0.22, y=1.12)
    lv = ["L3", "L6"]
    cats = [("lineages_with_>=2_annotated_plasmids", "≥ 2 annotated\nplasmids"),
            ("spanning_>=2_sectors", "in ≥ 2\nsectors"),
            ("same_key_AMR_gene_in_>=2_sectors", "same key gene\nin ≥ 2 sectors")]
    x = np.arange(3); w = 0.36
    for k, (l, col) in enumerate(zip(lv, [RAMP[3], RAMP[6]])):
        vals = [oh[l][c] for c, _ in cats]
        ax.bar(x + (k - 0.5) * w, vals, w - 0.04, color=col, label=f"{l} clusters", zorder=3)
        for xi, v in zip(x + (k - 0.5) * w, vals):
            ax.text(xi, v * 1.08, f"{v:,}", ha="center", va="bottom", fontsize=5.5, color=INK2)
    ax.set_yscale("log"); ax.set_ylim(50, 30000)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_xticks(x); ax.set_xticklabels([n for _, n in cats], fontsize=6)
    ax.set_ylabel("pLIN clusters (log)")
    ax.legend(frameon=False, loc="upper right", ncol=2)
    ax.set_title(f"{oh['matched_to_PLSDB_2024_05_31']:,} plasmids matched\nto PLSDB metadata",
                 loc="left", fontweight="bold")
    style(ax)

    # b: largest outbreak clones shared across sectors
    e = ex[ex.level == "L6"].sort_values("plasmids_with_gene", ascending=False).head(10).reset_index(drop=True)
    ax = fig.add_axes([0.585, 0.16, 0.15, 0.62]); label(ax, "b", x=-1.75, y=1.18)
    ys = np.arange(len(e))[::-1]
    for y, (_, r) in zip(ys, e.iterrows()):
        s = set(r.sectors.split(","))
        for j, sec in enumerate(SECTORS):
            ax.scatter(j, y, s=30, color=PLIN if sec in s else "#ffffff", edgecolor=PLIN if sec in s else GRID,
                       lw=0.8, zorder=3)
    ax.set_xticks(range(len(SECTORS))); ax.set_xticklabels(SECTORS, rotation=90)
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
    ax.set_yticks(ys); ax.set_yticklabels([f"{r.gene}  {r.lineage}" for _, r in e.iterrows()], fontsize=6)
    ax.set_xlim(-0.6, len(SECTORS) - 0.4); ax.set_ylim(-0.6, len(e) - 0.4)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    fig.text(0.36, 1.0, "Largest L6 outbreak clones found in \u2265 2 sectors (key resistance gene, pLIN code)",
             fontweight="bold", va="top")

    ax2 = fig.add_axes([0.76, 0.16, 0.10, 0.62]); ax2.set_ylim(-0.6, len(e) - 0.4)
    ax2.barh(ys, e.plasmids_with_gene, 0.6, color=RAMP[4], zorder=3)
    for y, v in zip(ys, e.plasmids_with_gene):
        ax2.text(v + 8, y, str(v), va="center", fontsize=6)
    ax2.set_yticks([]); ax2.set_xlabel("Plasmids with\nthe gene"); ax2.set_xlim(0, e.plasmids_with_gene.max() * 1.3)
    style(ax2, grid="x")
    ax3 = fig.add_axes([0.88, 0.16, 0.10, 0.62]); ax3.set_ylim(-0.6, len(e) - 0.4)
    ax3.barh(ys, e.countries, 0.6, color=RAMP[2], zorder=3)
    for y, v in zip(ys, e.countries):
        ax3.text(v + 1.5, y, str(v), va="center", fontsize=6)
    ax3.set_yticks([]); ax3.set_xlabel("Countries"); ax3.set_xlim(0, e.countries.max() * 1.35)
    style(ax3, grid="x")
    save(fig, "Fig5")


# ----------------------------------------------------------------------------- Supplementary
def sfig1():
    bs = js("build_stats.json")
    fig = plt.figure(figsize=(180 * MM, 120 * MM))

    # a: sources
    ax = fig.add_axes([0.08, 0.66, 0.88, 0.22]); label(ax, "a", x=-0.07, y=1.15)
    parts = [("PLSDB 2024_05_31", bs["source_PLSDB_2024_05_31"], RAMP[5]),
             ("replicon training set, not in PLSDB", bs["source_training_not_in_PLSDB"], RAMP[3]),
             ("other NCBI plasmids", bs["source_other_NCBI"], RAMP[1])]
    left = 0
    for nm, v, c in parts:
        ax.barh(0, v, 0.5, left=left, color=c, edgecolor=SURF, lw=1.5)
        if v > 10000:
            ax.text(left + v / 2, 0, f"{v:,}", ha="center", va="center", fontsize=6.5, color="#ffffff" if c != RAMP[1] else INK)
            ax.text(left + v / 2, 0.42, nm, ha="center", va="bottom", fontsize=6.5)
        else:
            ax.annotate(f"{nm}: {v:,}", (left + v / 2, 0.25), xytext=(left + v / 2, 0.62), ha="center", va="bottom",
                        fontsize=6.5, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
        left += v
    ax.set_xlim(0, left); ax.set_ylim(-0.4, 0.9); ax.set_yticks([])
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_xlabel("Unique plasmids in release")
    ax.set_title(f"Database db-2026.10.03: {bs['unique_accessions']:,} unique plasmids "
                 f"({bs['accessions_listed_twice']:,} accessions listed twice in the previous table merged)",
                 loc="left", fontweight="bold", pad=18)
    style(ax, grid=None); ax.spines["left"].set_visible(False)

    # b: processing
    ax = fig.add_axes([0.0, 0.08, 0.50, 0.46]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    label(ax, "b", x=0.0, y=1.0)
    steps = [(f"{bs['unique_accessions']:,} plasmid sequences", None),
             (f"{bs['proteins_predicted_all_ids']:,} predicted proteins\n(pyrodigal, metagenomic mode)",
              f"{bs['plasmids_without_proteins']:,} plasmids without a protein"),
             (f"{bs['unique_protein_sequences']:,} unique protein sequences", None),
             (f"{bs['protein_families_in_use']:,} protein families\n(MMseqs2, ≥ 50% identity, ≥ 80% coverage;\nidentical sequences in one family)",
              f"median {bs['protein_families_per_plasmid_median']:.0f} families per plasmid"),
             (f"{bs['kmer_hashes_total']:,} k-mer hashes in FracMinHash sketches\n(k = 21)",
              f"median {bs['kmer_hashes_per_plasmid_median']} per plasmid")]
    for k, (t, note) in enumerate(steps):
        y = 0.86 - k * 0.19
        box(ax, 0.06, y - 0.07, 0.62, 0.14, t, fs=6.3)
        if note:
            ax.text(0.70, y, note, fontsize=6, color=INK2, va="center")
        if k < len(steps) - 1:
            arrow(ax, 0.37, y - 0.07, 0.37, y - 0.12)

    # c: clusters per level
    ax = fig.add_axes([0.62, 0.10, 0.35, 0.40]); label(ax, "c", x=-0.24, y=1.1)
    cl = bs["clusters_per_level"]; lv = list(cl)
    ax.bar(range(6), [cl[l] for l in lv], 0.62, color=LEVEL_RAMP, zorder=3)
    for k, l in enumerate(lv):
        ax.text(k, cl[l] + 1500, f"{cl[l]:,}", ha="center", va="bottom", fontsize=5.8)
    ax.set_xticks(range(6)); ax.set_xticklabels(lv); ax.set_ylabel("Clusters")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_ylim(0, max(cl.values()) * 1.15)
    ax.set_title("Clusters per level", loc="left", fontweight="bold")
    style(ax)
    save(fig, "SFig1")


def sfig2():
    fig = plt.figure(figsize=(180 * MM, 130 * MM))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.02, 0.97, "Assignment of a new plasmid (identical in database build, app and command line)",
            fontsize=7.5, fontweight="bold", va="top")
    # inputs
    box(ax, 0.02, 0.74, 0.27, 0.15, "Protein families\npyrodigal genes; family of the identical\ncatalogue protein, else nearest\n"
        "(MMseqs2 ≥ 50% id, ≥ 80% cov), else new", fs=6.2)
    box(ax, 0.02, 0.54, 0.27, 0.15, "k-mer sketch\ncanonical 21-mers, splitmix64 hash;\nFracMinHash scale = power of two\n"
        "keeping ≥ 400 k-mers (max 256)", fs=6.2)
    box(ax, 0.02, 0.27, 0.27, 0.22, "Similarities to database plasmids\nprot: shared families / families of\nthe smaller "
        "plasmid (≥ min(3, size) shared)\nkcont: k-mers of the smaller plasmid\nfound in the other\n"
        "kmin: shared k-mers / k-mers of the\nlarger plasmid (symmetric)", fs=6.2, ha="left")
    # decision
    box(ax, 0.36, 0.60, 0.25, 0.16, "Most similar database plasmid\n(kmin; ties: earliest)\nkmin ≥ 0.80 ?",
        fc="#e8f1fc", ec=PLIN)
    box(ax, 0.70, 0.74, 0.28, 0.17, "Lineage path\ncopy its L1 to L5;\ninside that L5, nearest plasmid with\n"
        "kmin ≥ 0.95 gives L6, else a new L6", fc="#e8f1fc", ec=PLIN, fs=6.2)
    box(ax, 0.70, 0.30, 0.28, 0.34, "Founder path (L1 to L4)\nL1: earliest founder with prot ≥ 0.40\n"
        "L2: inside L1, earliest with prot ≥ 0.60\nL3: inside L2, earliest with kcont ≥ 0.50\n"
        "L4: inside L3, earliest with kcont ≥ 0.80\nno founder reached: the plasmid becomes\na new founder (new ID) at that level\n"
        "and all finer levels; L5 and L6 new\nno proteins: L1 = L2 = 0", fc="#e8f1fc", ec=PLIN, fs=6.2, ha="left")
    arrow(ax, 0.29, 0.81, 0.36, 0.70); arrow(ax, 0.29, 0.61, 0.36, 0.66); arrow(ax, 0.29, 0.38, 0.36, 0.62)
    arrow(ax, 0.61, 0.71, 0.70, 0.81); arrow(ax, 0.61, 0.64, 0.70, 0.50)
    ax.text(0.655, 0.79, "yes", fontsize=6, color=INK2, ha="center"); ax.text(0.655, 0.55, "no", fontsize=6, color=INK2, ha="center")
    # release guarantees
    box(ax, 0.02, 0.03, 0.96, 0.18,
        "Why codes are permanent\n"
        "• The release database is built once, in accession order; each database plasmid keeps its code in every later release.\n"
        "• New plasmids are added after all database plasmids, so they can join existing clusters or found new ones but never move an existing one.\n"
        "• Founders are fixed: a cluster is defined by its first member, not by a centroid or by single linkage, so clusters never merge.\n"
        "• A sequence already in the database is recognised by its exact hash and receives its published code.",
        fs=6.2, ha="left", fc="#f4f3f0")
    save(fig, "SFig2")


def sfig4():
    tm = tsv("confirm", "test_metrics.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(180 * MM, 95 * MM), sharey=True)
    order = [m for m in tm.method.unique()]
    order = sorted(order, key=lambda m: (0 if m.startswith("pLIN v4.1") else 1 if m.startswith("pLIN v4 ") else
                                         2 if m.startswith("MOB") else 3 if m.startswith("pling") else 4, m))
    col = lambda m: PLIN if m.startswith("pLIN v4.1") else LEGACY if m.startswith("pLIN v4 ") else \
        MOB if m.startswith("MOB") else PLING if m.startswith("pling") else MGE
    ys = np.arange(len(order))[::-1]
    for ax, truth, ttl, lt in zip(axes, ["same_lineage", "related_backbone"], ["Same lineage", "Related backbone"], "ab"):
        for y, m in zip(ys, order):
            r = metric(tm, m, truth, None)
            ax.plot([r.F1_w_lo, r.F1_w_hi], [y, y], color=col(m), lw=1.2, zorder=3)
            ax.plot(r.F1_w, y, "o", color=col(m), ms=4.5, mec=SURF, mew=0.8, zorder=4)
            ax.text(1.02, y, f"{r.F1_w:.2f}  P {r.precision_w:.2f}  R {r.recall_w:.2f}", va="center", fontsize=5.6,
                    color=INK2, transform=ax.get_yaxis_transform())
        ax.set_xlim(0, 1); ax.set_title(ttl, loc="left", fontweight="bold"); ax.set_xlabel("Weighted F1 (95% CI)")
        style(ax, grid="x"); label(ax, lt, x=-0.02 if lt == "b" else -0.35)
    axes[0].set_yticks(ys); axes[0].set_yticklabels(order)
    fig.subplots_adjust(wspace=0.55, left=0.17, right=0.86)
    save(fig, "SFig4")


def sfig5():
    oc, loo = tsv("case_studies", "outbreak_comparators.tsv"), tsv("case_studies", "loo_metrics.tsv")
    om = tsv("case_studies", "outbreak_metrics.tsv")
    fig = plt.figure(figsize=(180 * MM, 90 * MM))
    ax = fig.add_axes([0.20, 0.14, 0.30, 0.74]); label(ax, "a", x=-0.55, y=1.04)
    rows = [f"pLIN v4.1 L{k}" for k in range(1, 7)] + [m for m in oc.method if not m.startswith("pLIN")]
    allm = pd.concat([oc[["method", "same_study_together_pct", "different_study_together_pct"]],
                      om.assign(method="pLIN v4.1 " + om.level).rename(columns={})[
                          ["method", "same_study_together_pct", "different_study_together_pct"]]]).drop_duplicates("method")
    allm = allm.set_index("method").loc[rows]
    ys = np.arange(len(rows))[::-1]
    c = [PLIN if m.startswith("pLIN") else MOB if m.startswith("MOB") else PLING if m.startswith("pling") else MGE for m in rows]
    ax.barh(ys + 0.19, allm.same_study_together_pct, 0.36, color=c, zorder=3)
    ax.barh(ys - 0.19, allm.different_study_together_pct, 0.36, color=c, alpha=0.45, zorder=3)
    for y, a, b in zip(ys, allm.same_study_together_pct, allm.different_study_together_pct):
        ax.text(a + 1, y + 0.19, f"{a:.1f}", va="center", fontsize=5.5)
        ax.text(b + 1, y - 0.19, f"{b:.2f}", va="center", fontsize=5.5, color=INK2)
    ax.set_yticks(ys); ax.set_yticklabels(rows); ax.set_xlim(0, 118)
    ax.set_xlabel("Pairs grouped together (%)\nsolid: same outbreak; light: different outbreaks")
    ax.set_title("Comparative outbreak study (all levels)", loc="left", fontweight="bold")
    style(ax, grid="x")

    ax = fig.add_axes([0.62, 0.14, 0.36, 0.74]); label(ax, "b", x=-0.17, y=1.04)
    lv = np.arange(1, 7)
    series = [("recovery_pct", RAMP[3], "own code recovered"), ("same_study_linked_pct", RAMP[6], "linked to its own outbreak"),
              ("false_link_pct", MUTED, "linked to another outbreak")]
    for colname, cc, nm in series:
        ax.plot(lv, loo[colname], "-o", color=cc, lw=1.4, ms=4, mec=SURF, mew=0.8, zorder=3)
        ax.text(6.15, loo[colname].iloc[-1], nm, fontsize=6, va="center")
    ax.set_xticks(lv); ax.set_xticklabels([f"L{i}" for i in lv]); ax.set_xlim(0.7, 8.6); ax.set_ylim(0, 105)
    ax.set_ylabel("% of plasmids")
    ax.set_title(f"Leave-one-out ({int(loo.plasmids_in_database.iloc[0])} plasmids; "
                 f"{int(loo.plasmids_with_outbreak_partner.iloc[0])} with an outbreak partner)", loc="left", fontweight="bold")
    style(ax)
    save(fig, "SFig5")


def sfig6():
    arr = tsv("case_studies", "prospective_arrivals.tsv")
    arr["date"] = pd.to_datetime(arr.date)
    fig, axes = plt.subplots(1, 2, figsize=(180 * MM, 95 * MM), sharey=True)
    st = arr.study.tolist()
    for ax, (col, nm, c, lt) in zip(axes, [("pLIN_v41_at_arrival", "pLIN v4.1 L5", PLIN, "a"),
                                          ("MOB_secondary", "MOB-suite secondary", MOB, "b")]):
        lab = [".".join(x.split(".")[:5]) for x in arr[col]] if col.startswith("pLIN") else \
            [None if pd.isna(x) else x for x in arr[col]]
        for i in range(len(arr)):
            same = [j for j in range(i) if st[j] == st[i]]
            other = [j for j in range(i) if st[j] != st[i]]
            if not same:
                ax.scatter(arr.date[i], i + 1, s=9, color=GRID, zorder=2)
                continue
            warn = lab[i] is not None and any(lab[j] == lab[i] for j in same)
            alarm = lab[i] is not None and any(lab[j] == lab[i] for j in other)
            ax.scatter(arr.date[i], i + 1, s=18, color=c if warn else "#ffffff", edgecolor=c, lw=0.9, zorder=3)
            if alarm:
                ax.scatter(arr.date[i], i + 1, s=34, marker="x", color=INK, lw=0.8, zorder=4)
        ax.set_title(nm, loc="left", fontweight="bold"); ax.set_xlabel("NCBI deposition date")
        style(ax, grid="x"); label(ax, lt, x=-0.12 if lt == "a" else -0.04)
    axes[0].set_ylabel("Arrival order")
    axes[1].text(1.0, -0.2, "filled: early warning (linked to an earlier isolate of its outbreak); open: missed;\n"
                 "x: false alarm (linked to an earlier isolate of another outbreak); grey: first isolate of its outbreak",
                 transform=axes[1].transAxes, ha="right", va="top", fontsize=6, color=INK2)
    save(fig, "SFig6")


def sfig7():
    p = tsv("ptu_agreement.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(180 * MM, 70 * MM), sharey=True)
    for ax, (s, lt) in zip(axes, [("all release plasmids with a PTU", "a"), ("confirmatory test set", "b")]):
        d = p[p.set == s]
        lv = d[d.method.str.startswith("pLIN")]
        ax.plot(range(1, 7), lv.ARI, "-o", color=PLIN, lw=1.5, ms=4, mec=SURF, mew=0.8, zorder=3, label="pLIN v4.1 levels")
        for k, (m, ls) in enumerate([("MOB-suite primary", "-"), ("MOB-suite secondary", (0, (3, 2)))]):
            v = d[d.method == m].ARI.iloc[0]
            ax.axhline(v, color=MOB, lw=1.0, ls=ls, zorder=2)
            ax.text(8.35, v + 0.012, f"{m} {v:.2f}", fontsize=5.8, va="bottom", ha="right")
        ax.annotate(f"{lv.ARI.iloc[0]:.2f}", (1, lv.ARI.iloc[0]), xytext=(0, 6), textcoords="offset points", ha="center", fontsize=6)
        ax.set_xticks(range(1, 7)); ax.set_xticklabels([f"L{i}" for i in range(1, 7)]); ax.set_xlim(0.7, 8.4)
        ax.set_ylim(0, 1); ax.set_title(f"{s[0].upper() + s[1:]} (n = {int(d.n.iloc[0]):,})", loc="left", fontweight="bold")
        style(ax); label(ax, lt, x=-0.15 if lt == "a" else -0.05)
    axes[0].set_ylabel("Adjusted Rand index with PTU")
    save(fig, "SFig7")


def graphical_abstract():
    """Graphical abstract (NAR Genomics and Bioinformatics): the code levels and the headline results."""
    fl = json.load(open(os.path.join(BASE_DIR, "docs", "PLIN_FACTS.json")))
    code = fl["oh_oxa48_code"]["text"].split(".")
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.04, 0.92, "pLIN: one permanent code per plasmid, from backbone family to outbreak clone",
            fontsize=10, fontweight="bold", va="top")
    w, g = 0.143, 0.012
    for k, (lv, meaning, *_ ) in enumerate(LEVEL_TABLE):
        x = 0.04 + k * (w + g)
        ax.add_patch(FancyBboxPatch((x, 0.58), w, 0.16, boxstyle="round,pad=0,rounding_size=0.02", fc=LEVEL_RAMP[k], ec="none"))
        ax.text(x + w / 2, 0.66, code[k], ha="center", va="center", fontsize=12, fontweight="bold",
                color="#ffffff" if k else INK)
        ax.text(x + w / 2, 0.52, f"{lv} {meaning}", ha="center", va="center", fontsize=7)
    ax.text(0.04, 0.44, "protein families and k-mers, fixed founders", fontsize=7.5, color=INK2)
    ax.text(0.96, 0.44, "nearest-relative lineage", fontsize=7.5, color=INK2, ha="right")
    stats = [(fl["lineage_F1"]["text"], f"lineage F1 (pre-registered test)\nMOB-suite {fl['lineage_F1_MOB']['text']}"),
             (fl["stability_retained"]["text"], "codes unchanged as\nthe database grows"),
             (f"{fl['requery_n']['text']}/{fl['requery_n']['text']}", "codes reproduced\nfrom sequence"),
             (fl["onehealth_clones_multisector"]["text"], "plasmid clones shared\nbetween One Health sectors")]
    for k, (big, small) in enumerate(stats):
        x = 0.04 + k * 0.24
        ax.text(x, 0.27, big, fontsize=18, fontweight="bold", color=PLIN, va="center")
        ax.text(x, 0.11, small, fontsize=7.5, color=INK2, va="center")
    fig.savefig(os.path.join(OUT, "Graphical_abstract.png"), dpi=300)
    plt.close(fig)
    print("wrote Graphical_abstract")


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, fig5, sfig1, sfig2, sfig4, sfig5, sfig6, sfig7, graphical_abstract):
        f()
