#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Sensitivity of the confirmatory results to the truth definitions (secondary analysis, not pre-registered).

The registered truth: same lineage = AF_min >= 0.80 and ANI >= 99%; related backbone = backbone
AF of the shorter plasmid >= 0.50. Because pLIN's lineage threshold is also 0.80 (symmetric k-mer
similarity), we re-score the same 3,337 confirmatory pairs (same weights, same labels, same bootstrap)
under other definitions:
  same lineage      AF_min in {0.6, 0.7, 0.8, 0.9} x ANI in {98, 99, 99.5}
  related backbone  backbone AF in {0.3, 0.4, 0.5, 0.6, 0.7}
For each definition: weighted F1 of every method at every level, the best level of each method, and the
registered comparisons (pLIN L5 vs MOB-suite secondary; pLIN L1 vs MOB-suite primary) with 95% bootstrap
CIs of the difference.

Usage:
  python v41_truth_sensitivity.py
Output: output/backbone_v41/confirm/truth_sensitivity.tsv, truth_sensitivity_best.tsv
"""

import os

import numpy as np
import pandas as pd

from v41_confirm_eval import C, labels
from v4_evaluate import N_BOOT, SEED, boot_f1, metrics, same_group


def main():
    pairs = pd.read_csv(os.path.join(C, "truth_pairs.tsv"), sep="\t").reset_index(drop=True)
    a, b, w = pairs.plasmid_shorter, pairs.plasmid_longer, pairs.weight.values
    ids = sorted(set(a) | set(b))
    lab = labels(sorted(pd.read_csv(os.path.join(C, "test_plasmids.tsv"), sep="\t").plasmid_id))
    pos = {i: k for k, i in enumerate(ids)}
    rng = np.random.default_rng(SEED)
    counts = rng.multinomial(len(ids), np.full(len(ids), 1 / len(ids)), size=N_BOOT)
    mult = counts[:, [pos[x] for x in a]] * counts[:, [pos[x] for x in b]]
    preds = {name: same_group(l, a, b) for name, l in lab.items()}

    af, ani, bb = pairs.AF_min.fillna(0).values, pairs.ANI_aln.fillna(0).values, pairs.bb_AF_shorter.fillna(0).values
    truths = []
    for t_af in (0.6, 0.7, 0.8, 0.9):
        for t_ani in (98.0, 99.0, 99.5):
            truths.append(("same_lineage", f"AF>={t_af}, ANI>={t_ani}", (af >= t_af) & (ani >= t_ani),
                           t_af == 0.8 and t_ani == 99.0))
    for t_bb in (0.3, 0.4, 0.5, 0.6, 0.7):
        truths.append(("related_backbone", f"backbone AF>={t_bb}", bb >= t_bb, t_bb == 0.5))

    # the registered definitions must reproduce the registered truth exactly
    for kind, name, y, reg in truths:
        if reg:
            assert (y == pairs[kind].values).all(), f"registered {kind} truth not reproduced"

    rows, best = [], []
    for kind, name, y, reg in truths:
        boots = {}
        for m, p in preds.items():
            pr, rc, f1 = metrics(p, y, w)
            boots[m] = boot_f1(p, y, w, mult)
            rows.append({"truth": kind, "definition": name, "registered": reg, "positive_pairs": int(y.sum()),
                         "method": m, "precision_w": pr, "recall_w": rc, "F1_w": f1})
        res = pd.DataFrame([r for r in rows if r["definition"] == name])
        tool = res.method.str.replace(r" (L\d|primary|secondary|community|subcommunity)$", "", regex=True)
        for t, g in res.groupby(tool):
            top = g.sort_values("F1_w", ascending=False).iloc[0]
            best.append({"truth": kind, "definition": name, "registered": reg, "tool": t, "best_level": top.method,
                         "F1_w": top.F1_w})
        ref = ("pLIN v4.1 L5", "MOB-suite secondary") if kind == "same_lineage" else ("pLIN v4.1 L1", "MOB-suite primary")
        d = boots[ref[0]] - boots[ref[1]]
        f = {m: float(res[res.method == m].F1_w.iloc[0]) for m in ref}
        best.append({"truth": kind, "definition": name, "registered": reg, "tool": "registered comparison",
                     "best_level": f"{ref[0]} minus {ref[1]}", "F1_w": f[ref[0]] - f[ref[1]],
                     "diff_lo": float(np.percentile(d, 2.5)), "diff_hi": float(np.percentile(d, 97.5))})
    pd.DataFrame(rows).to_csv(os.path.join(C, "truth_sensitivity.tsv"), sep="\t", index=False)
    bt = pd.DataFrame(best)
    bt.to_csv(os.path.join(C, "truth_sensitivity_best.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 250)
    piv = bt[bt.tool != "registered comparison"].pivot_table(index=["truth", "definition"], columns="tool", values="F1_w")
    print(piv.round(2).to_string())
    print(bt[bt.tool == "registered comparison"][["truth", "definition", "F1_w", "diff_lo", "diff_hi"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
