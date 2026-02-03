#!/usr/bin/env python3
"""
pLIN Test Cladogram
Generates a publication-quality cladogram (dendrogram) for test plasmids
with pLIN code annotations and threshold markers.
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
from itertools import product as iter_product
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from Bio import SeqIO

# ── Configuration ──────────────────────────────────────────────────────────────

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.join(BASE_DIR, "test_plasmids")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "test")

PLIN_THRESHOLDS = {
    "A": 0.150,
    "B": 0.100,
    "C": 0.050,
    "D": 0.020,
    "E": 0.010,
    "F": 0.001,
}

THRESHOLD_COLORS = {
    "A": "#E53935",
    "B": "#FB8C00",
    "C": "#FDD835",
    "D": "#43A047",
    "E": "#1E88E5",
    "F": "#8E24AA",
}

PLIN_LEVEL_NAMES = {
    "A": "Family",
    "B": "Subfamily",
    "C": "Cluster",
    "D": "Subcluster",
    "E": "Clone",
    "F": "Strain",
}

STRAIN_COLORS = [
    "#2196F3", "#FF9800", "#4CAF50", "#E91E63", "#9C27B0",
    "#00BCD4", "#FF5722", "#795548", "#607D8B", "#CDDC39",
    "#F44336", "#3F51B5", "#009688", "#FFC107", "#8BC34A",
]


# ── Data loading & computation ────────────────────────────────────────────────

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


# ── Figure 1: Rectangular cladogram ──────────────────────────────────────────

def plot_rectangular_cladogram(Z, labels, plin_codes, strain_clusters):
    strain_cmap = get_strain_cmap(strain_clusters)
    unique_strains = sorted(set(strain_clusters))

    fig, (ax_dendro, ax_plin) = plt.subplots(
        1, 2, figsize=(14, 10),
        gridspec_kw={"width_ratios": [3, 2], "wspace": 0.02},
    )

    # -- Dendrogram panel --
    with plt.rc_context({"lines.linewidth": 2.0}):
        ddata = dendrogram(
            Z, labels=labels, orientation="right",
            leaf_font_size=9.5, ax=ax_dendro,
            color_threshold=0, above_threshold_color="#444444",
        )

    # Clamp x-axis to data range (max merge distance + 10% padding)
    max_merge = Z[:, 2].max()
    x_max = max_merge * 1.15
    ax_dendro.set_xlim(0, x_max)

    # Colour leaf labels
    leaf_colors = {labels[i]: strain_cmap[sc] for i, sc in enumerate(strain_clusters)}
    for lbl in ax_dendro.get_yticklabels():
        txt = lbl.get_text()
        if txt in leaf_colors:
            lbl.set_color(leaf_colors[txt])
            lbl.set_fontweight("bold")
            lbl.set_fontsize(10)

    # Threshold lines (only those within visible range)
    for bname, thresh in PLIN_THRESHOLDS.items():
        if thresh <= x_max:
            color = THRESHOLD_COLORS[bname]
            ax_dendro.axvline(x=thresh, color=color, linestyle="--",
                              linewidth=1.3, alpha=0.7, zorder=0)
            ax_dendro.text(
                thresh, ax_dendro.get_ylim()[1] * 1.01,
                f" {bname} ({PLIN_LEVEL_NAMES[bname]}) d\u2264{thresh}",
                fontsize=7.5, color=color, fontweight="bold",
                ha="left", va="bottom", rotation=0,
            )

    # Add threshold legend for those beyond visible range
    offscreen = [(b, t) for b, t in PLIN_THRESHOLDS.items() if t > x_max]
    if offscreen:
        text_lines = ["Thresholds beyond range:"]
        for b, t in offscreen:
            text_lines.append(f"  {b} ({PLIN_LEVEL_NAMES[b]}): d\u2264{t}")
        ax_dendro.text(
            0.98, 0.02, "\n".join(text_lines),
            transform=ax_dendro.transAxes,
            fontsize=7, va="bottom", ha="right",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.9),
        )

    ax_dendro.set_xlabel("Cosine Distance (4-mer composition)", fontsize=11, fontweight="bold")
    ax_dendro.set_title("Dendrogram", fontsize=12, fontweight="bold")

    # -- pLIN annotation panel --
    ax_plin.set_ylim(ax_dendro.get_ylim())
    ax_plin.set_xlim(0, 1)
    ax_plin.axis("off")
    ax_plin.set_title("pLIN Code", fontsize=12, fontweight="bold")

    leaf_order = ddata["leaves"]
    y_positions = [5 + i * 10 for i in range(len(leaf_order))]
    label_to_plin = dict(zip(labels, plin_codes))
    label_to_strain = dict(zip(labels, strain_clusters))

    for i, leaf_idx in enumerate(leaf_order):
        lbl = labels[leaf_idx]
        plin = label_to_plin[lbl]
        sc = label_to_strain[lbl]
        color = strain_cmap[sc]
        y = y_positions[i]
        ax_plin.text(0.05, y, plin, fontsize=9.5, fontfamily="monospace",
                     fontweight="bold", color=color, va="center")

    # Strain legend
    legend_elements = [
        mpatches.Patch(color=strain_cmap[sc],
                       label=f"Strain {sc} (n={sum(1 for s in strain_clusters if s == sc)})")
        for sc in unique_strains
    ]
    ax_dendro.legend(
        handles=legend_elements, title="pLIN Strain (F)",
        loc="upper right", fontsize=8, title_fontsize=9, framealpha=0.9,
    )

    fig.suptitle(
        "pLIN Cladogram \u2014 IncX Test Plasmids (n=22)\n"
        "Single-linkage clustering on tetranucleotide composition",
        fontsize=14, fontweight="bold", y=0.98,
    )

    return fig


# ── Figure 2: Circular cladogram ────────────────────────────────────────────

def plot_circular_cladogram(Z, labels, plin_codes, strain_clusters):
    strain_cmap = get_strain_cmap(strain_clusters)
    unique_strains = sorted(set(strain_clusters))

    fig, ax = plt.subplots(figsize=(14, 14), subplot_kw={"polar": True})

    n = len(labels)
    ddata = dendrogram(Z, labels=labels, no_plot=True)
    leaf_order = ddata["leaves"]

    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    angle_map = {leaf_order[i]: angles[i] for i in range(n)}

    max_dist = max(Z[:, 2]) * 1.15
    radius_scale = 0.85

    node_positions = {}
    for i in range(n):
        node_positions[i] = (angle_map[i], 0.0)

    for idx, (c1, c2, dist, count) in enumerate(Z):
        c1, c2 = int(c1), int(c2)
        a1, d1 = node_positions[c1]
        a2, d2 = node_positions[c2]

        node_id = n + idx
        avg_angle = np.arctan2(
            (np.sin(a1) + np.sin(a2)) / 2,
            (np.cos(a1) + np.cos(a2)) / 2,
        )
        if avg_angle < 0:
            avg_angle += 2 * np.pi
        node_positions[node_id] = (avg_angle, dist)

        r1_inner = (1.0 - d1 / max_dist) * radius_scale
        r1_outer = (1.0 - dist / max_dist) * radius_scale
        r2_inner = (1.0 - d2 / max_dist) * radius_scale
        r2_outer = (1.0 - dist / max_dist) * radius_scale

        ax.plot([a1, a1], [r1_inner, r1_outer], color="#555555", linewidth=1.5)
        ax.plot([a2, a2], [r2_inner, r2_outer], color="#555555", linewidth=1.5)

        r_merge = (1.0 - dist / max_dist) * radius_scale
        a_min, a_max = min(a1, a2), max(a1, a2)
        if (a_max - a_min) > np.pi:
            arc_angles = np.linspace(a_max, a_min + 2 * np.pi, 50)
        else:
            arc_angles = np.linspace(a_min, a_max, 50)
        ax.plot(arc_angles, np.full_like(arc_angles, r_merge), color="#555555", linewidth=1.5)

    # Leaf markers and labels
    for i in range(n):
        angle = angle_map[i]
        r = radius_scale + 0.02
        sc = strain_clusters[i]
        color = strain_cmap[sc]

        ax.scatter(angle, r, s=80, c=color, zorder=5, edgecolors="white", linewidth=0.5)

        label_r = radius_scale + 0.07
        rotation = np.degrees(angle) - 90
        ha = "left"
        if 90 < np.degrees(angle) < 270:
            rotation += 180
            ha = "right"

        ax.text(
            angle, label_r,
            f"{labels[i]}  [{plin_codes[i]}]",
            fontsize=7.5, fontfamily="monospace",
            fontweight="bold", color=color,
            ha=ha, va="center",
            rotation=rotation, rotation_mode="anchor",
        )

    # Threshold circles (only E and F are within data range typically)
    for bname, thresh in PLIN_THRESHOLDS.items():
        r_thresh = (1.0 - thresh / max_dist) * radius_scale
        if 0.05 < r_thresh < radius_scale:
            theta_circle = np.linspace(0, 2 * np.pi, 200)
            ax.plot(theta_circle, np.full(200, r_thresh),
                    color=THRESHOLD_COLORS[bname], linestyle="--",
                    linewidth=1.0, alpha=0.5, label=f"{bname} d\u2264{thresh}")

    ax.set_ylim(0, radius_scale + 0.22)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines["polar"].set_visible(False)
    ax.grid(False)

    ax.set_title(
        "pLIN Circular Cladogram \u2014 IncX Test Plasmids (n=22)\n"
        "Tetranucleotide composition | Single-linkage clustering",
        fontsize=14, fontweight="bold", pad=30,
    )

    legend_elements = [
        mpatches.Patch(color=strain_cmap[sc], label=f"Strain {sc}")
        for sc in unique_strains
    ]
    ax.legend(
        handles=legend_elements, title="pLIN Strain (F)",
        loc="lower left", bbox_to_anchor=(-0.05, -0.05),
        fontsize=8, title_fontsize=9, framealpha=0.9,
    )

    return fig


# ── Figure 3: Cladogram + heatmap ───────────────────────────────────────────

def plot_cladogram_with_heatmap(Z, labels, plin_codes, cluster_assignments, records):
    fig = plt.figure(figsize=(18, 10))

    # Layout: dendrogram | heatmap | metadata
    ax_dendro = fig.add_axes([0.02, 0.08, 0.25, 0.82])
    ax_heat = fig.add_axes([0.30, 0.08, 0.28, 0.82])
    ax_meta = fig.add_axes([0.62, 0.08, 0.35, 0.82])

    # -- Dendrogram --
    ddata = dendrogram(
        Z, labels=labels, orientation="right",
        leaf_font_size=8.5, ax=ax_dendro,
        color_threshold=0, above_threshold_color="#555555",
    )

    # Clamp x-axis
    max_merge = Z[:, 2].max()
    ax_dendro.set_xlim(0, max_merge * 1.15)

    # Only show threshold lines within visible range
    for bname, thresh in PLIN_THRESHOLDS.items():
        if thresh <= max_merge * 1.15:
            ax_dendro.axvline(x=thresh, color=THRESHOLD_COLORS[bname],
                              linestyle="--", linewidth=1.0, alpha=0.6)

    ax_dendro.set_xlabel("Cosine Distance", fontsize=10)
    ax_dendro.set_title("Dendrogram", fontsize=11, fontweight="bold")

    leaf_order = ddata["leaves"]

    # -- Heatmap --
    bin_labels = list(PLIN_THRESHOLDS.keys())
    heat_data = np.zeros((len(labels), len(bin_labels)))
    for j, bname in enumerate(bin_labels):
        for i, leaf_idx in enumerate(leaf_order):
            heat_data[i, j] = cluster_assignments[bname][leaf_idx]

    im = ax_heat.imshow(heat_data, aspect="auto", cmap="tab20", interpolation="nearest")

    ax_heat.set_xticks(range(len(bin_labels)))
    ax_heat.set_xticklabels([f"{b}\n({PLIN_LEVEL_NAMES[b]})" for b in bin_labels], fontsize=9)
    ax_heat.set_yticks(range(len(labels)))
    ax_heat.set_yticklabels([labels[i] for i in leaf_order], fontsize=8)

    for i in range(heat_data.shape[0]):
        for j in range(heat_data.shape[1]):
            val = int(heat_data[i, j])
            ax_heat.text(j, i, str(val), ha="center", va="center",
                         fontsize=8, fontweight="bold", color="white",
                         path_effects=[pe.withStroke(linewidth=2, foreground="black")])

    ax_heat.set_title("pLIN Cluster Assignments", fontsize=11, fontweight="bold")

    # -- Metadata panel --
    ax_meta.axis("off")

    col_headers = ["Plasmid", "Length (bp)", "pLIN Code"]
    header_y = len(labels)
    for j, h in enumerate(col_headers):
        ax_meta.text(j * 0.35, header_y + 0.3, h, fontsize=9,
                     fontweight="bold", ha="left", va="center")
    ax_meta.axhline(y=header_y, color="black", linewidth=0.8, xmin=0.0, xmax=0.95)

    for i, leaf_idx in enumerate(leaf_order):
        y = len(labels) - 1 - i
        rec = records[leaf_idx]
        plin = plin_codes[leaf_idx]
        ax_meta.text(0.0, y, rec["source_file"].replace(".fasta", ""),
                     fontsize=8, ha="left", va="center", fontfamily="monospace")
        ax_meta.text(0.35, y, f"{rec['length']:,}",
                     fontsize=8, ha="left", va="center")
        ax_meta.text(0.70, y, plin,
                     fontsize=8, ha="left", va="center",
                     fontfamily="monospace", fontweight="bold", color="#1565C0")

    ax_meta.set_xlim(-0.05, 1.1)
    ax_meta.set_ylim(-1, len(labels) + 1)
    ax_meta.set_title("Plasmid Metadata", fontsize=11, fontweight="bold")

    fig.suptitle(
        "pLIN Cladogram with Hierarchical Cluster Assignments \u2014 IncX Test Plasmids (n=22)",
        fontsize=14, fontweight="bold", y=0.97,
    )

    return fig


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("pLIN Cladogram Generator \u2014 Test Plasmids")
    print("=" * 70)

    print("\n[1/5] Loading sequences ...")
    records = load_sequences()
    print(f"  Loaded {len(records)} sequences")

    if len(records) < 2:
        print("  ERROR: Need at least 2 sequences.")
        return

    print("[2/5] Computing 4-mer vectors ...")
    vectors = compute_kmer_vectors(records, k=4)

    print("[3/5] Computing distances and clustering ...")
    dist_condensed = pdist(vectors, metric="cosine")
    Z = linkage(dist_condensed, method="single")

    plin_codes, cluster_assignments = get_plin_codes(records, Z)
    strain_clusters = list(cluster_assignments["F"])

    labels = [r["source_file"].replace(".fasta", "") for r in records]

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── Figure 1: Rectangular cladogram ──────────────────────────────────────
    print("[4/5] Generating rectangular cladogram ...")
    fig1 = plot_rectangular_cladogram(Z, labels, plin_codes, strain_clusters)
    for ext in ["png", "pdf"]:
        fig1.savefig(os.path.join(OUTPUT_DIR, f"test_cladogram_rectangular.{ext}"),
                     dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig1)
    print("  Saved: test_cladogram_rectangular.png/pdf")

    # ── Figure 2: Circular cladogram ─────────────────────────────────────────
    print("         Generating circular cladogram ...")
    fig2 = plot_circular_cladogram(Z, labels, plin_codes, strain_clusters)
    for ext in ["png", "pdf"]:
        fig2.savefig(os.path.join(OUTPUT_DIR, f"test_cladogram_circular.{ext}"),
                     dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig2)
    print("  Saved: test_cladogram_circular.png/pdf")

    # ── Figure 3: Cladogram + heatmap ────────────────────────────────────────
    print("[5/5] Generating cladogram with heatmap ...")
    fig3 = plot_cladogram_with_heatmap(Z, labels, plin_codes, cluster_assignments, records)
    for ext in ["png", "pdf"]:
        fig3.savefig(os.path.join(OUTPUT_DIR, f"test_cladogram_heatmap.{ext}"),
                     dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig3)
    print("  Saved: test_cladogram_heatmap.png/pdf")

    print(f"\n{'=' * 70}")
    print(f"All cladograms saved to: {OUTPUT_DIR}/")
    print("  - test_cladogram_rectangular.png/pdf")
    print("  - test_cladogram_circular.png/pdf")
    print("  - test_cladogram_heatmap.png/pdf")
    print("=" * 70)


if __name__ == "__main__":
    main()
