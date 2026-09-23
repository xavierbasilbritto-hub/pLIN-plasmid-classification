#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
"""
Compare pLIN's uniform, IncX-calibrated thresholds against per-Inc-group
adaptive thresholds derived independently from each group's own within-group
cosine-distance distribution.

This directly answers Reviewer 2's methodological objection: thresholds were
"calibrated on 178 IncX-like plasmids... and applied uniformly rather than
re-derived per group... Is this really a valid thing to do? IncX plasmids
are a bit atypical (pi protein replication etc)."

Uses plin_app.py::calibrate_inc_thresholds(), which already implements
per-group quantile-based adaptive calibration but was previously unused
outside the interactive Streamlit app (opt-in checkbox, off by default).

Output: a TSV summarising, for every one of the 28 Inc/Rep groups, how far
each of the 6 uniform (IncX-derived) thresholds sits from the equivalent
quantile-derived threshold for that group's own distance distribution.
"""

import os
import sys
import json
import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)

# Import the calibration function and constants directly from plin_app.py
# without triggering Streamlit page execution.
import importlib.util

spec = importlib.util.spec_from_file_location("plin_app_module", os.path.join(BASE_DIR, "plin_app.py"))
plin_app_module = importlib.util.module_from_spec(spec)

# plin_app.py runs Streamlit page-config / UI code at import time in some
# versions; guard against that by only importing if it exposes the function
# we need without side effects. We instead read the two relevant pieces of
# logic (thresholds + calibration) directly to avoid importing the whole
# Streamlit app.

CLASSIFIER_PATH = os.path.join(BASE_DIR, "data", "inc_classifier.npz")

PLIN_THRESHOLDS = {
    "A": 0.150,
    "B": 0.100,
    "C": 0.050,
    "D": 0.020,
    "E": 0.010,
    "F": 0.001,
}

LEVEL_QUANTILES = {
    "A": 0.99,
    "B": 0.95,
    "C": 0.75,
    "D": 0.50,
    "E": 0.25,
    "F": 0.05,
}

LEVEL_NAMES = {
    "A": "L1 (Family)",
    "B": "L2 (Subfamily)",
    "C": "L3 (Cluster)",
    "D": "L4 (Subcluster)",
    "E": "L5 (Clone group)",
    "F": "L6 (Lineage)",
}


def calibrate_inc_thresholds():
    """Per-Inc-group quantile-based threshold calibration.

    Reimplemented here (rather than imported from plin_app.py) to avoid
    importing the Streamlit application module for a standalone analysis
    script; logic is identical to plin_app.py::calibrate_inc_thresholds().
    """
    from scipy.spatial.distance import pdist

    data = np.load(CLASSIFIER_PATH, allow_pickle=True)
    X = data["X"]
    y = data["y"]
    group_names = [str(g) for g in data["group_names"]]

    calibrated = {}
    group_n = {}
    for i, name in enumerate(group_names):
        mask = y == i
        X_group = X[mask]
        group_n[name] = int(X_group.shape[0])
        if X_group.shape[0] < 10:
            continue

        n = X_group.shape[0]
        if n > 500:
            rng = np.random.default_rng(42)
            idx = rng.choice(n, 500, replace=False)
            X_sub = X_group[idx]
        else:
            X_sub = X_group

        dists = pdist(X_sub.astype(np.float64), metric="cosine")

        thresholds = {}
        for level, q in LEVEL_QUANTILES.items():
            val = float(np.quantile(dists, q))
            thresholds[level] = round(max(val, 0.0005), 6)

        levels = list("ABCDEF")
        for j in range(1, len(levels)):
            if thresholds[levels[j]] >= thresholds[levels[j - 1]]:
                thresholds[levels[j]] = round(thresholds[levels[j - 1]] * 0.7, 6)

        calibrated[name] = thresholds

    return calibrated, group_n


def main():
    print("Loading classifier training data and calibrating per-group thresholds ...")
    calibrated, group_n = calibrate_inc_thresholds()

    is_incx_like = {"IncX1", "IncX3", "IncX4"}

    rows = []
    for group, thresholds in calibrated.items():
        for level in "ABCDEF":
            uniform_val = PLIN_THRESHOLDS[level]
            adaptive_val = thresholds[level]
            pct_diff = (adaptive_val - uniform_val) / uniform_val * 100
            rows.append({
                "inc_group": group,
                "n_training": group_n.get(group, 0),
                "is_incx_like": group in is_incx_like,
                "level": level,
                "level_name": LEVEL_NAMES[level],
                "uniform_threshold": uniform_val,
                "adaptive_threshold": adaptive_val,
                "pct_difference": round(pct_diff, 1),
                "abs_pct_difference": round(abs(pct_diff), 1),
            })

    df = pd.DataFrame(rows)

    out_path = os.path.join(BASE_DIR, "output", "adaptive_threshold_comparison.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\nSaved: {out_path}  ({len(df)} rows, {df['inc_group'].nunique()} groups)")

    # Summary: how far off is the uniform (IncX-derived) threshold, on
    # average, for IncX-like groups themselves vs. all other groups?
    summary = (
        df.groupby(["is_incx_like", "level"])["abs_pct_difference"]
        .agg(["mean", "median", "max", "count"])
        .reset_index()
    )
    summary_path = os.path.join(BASE_DIR, "output", "adaptive_threshold_summary.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"Saved: {summary_path}")

    print("\n" + "=" * 78)
    print("SUMMARY: |% difference| between uniform (IncX-derived) and adaptive")
    print("(per-group quantile) thresholds, split by IncX-like vs other groups")
    print("=" * 78)
    print(summary.to_string(index=False))

    # Overall groups most divergent from the uniform scheme (by mean abs % diff
    # across all 6 levels) — this is the direct answer to "is applying IncX-
    # derived thresholds uniformly valid?"
    per_group = (
        df.groupby(["inc_group", "n_training", "is_incx_like"])["abs_pct_difference"]
        .mean()
        .reset_index()
        .sort_values("abs_pct_difference", ascending=False)
    )
    per_group_path = os.path.join(BASE_DIR, "output", "adaptive_threshold_per_group_divergence.tsv")
    per_group.to_csv(per_group_path, sep="\t", index=False)
    print(f"\nSaved: {per_group_path}")
    print("\nGroups ranked by mean |% divergence| from uniform IncX-derived thresholds:")
    print(per_group.to_string(index=False))

    with open(os.path.join(BASE_DIR, "output", "adaptive_threshold_comparison.meta.json"), "w") as fh:
        json.dump({
            "level_quantiles": LEVEL_QUANTILES,
            "uniform_thresholds": PLIN_THRESHOLDS,
            "incx_like_groups": sorted(is_incx_like),
            "n_groups_calibrated": len(calibrated),
        }, fh, indent=2)


if __name__ == "__main__":
    main()
