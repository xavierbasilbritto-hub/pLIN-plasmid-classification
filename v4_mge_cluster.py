#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
mge-cluster as a secondary comparator (pre-registration deviation 4).

mge-cluster 1.1.0, default settings (k = 31, min_cluster = 30, perplexity 30):
  model_all       --create on all 2,000 evaluation plasmids -> accuracy
  model_A         --create on snapshot A (seed-1 half)       -> labels before
  existing_A_all  --existing with model_A on all 2,000       -> update mode
(`--model_prefix mge-cluster` must be given explicitly in existing mode; the
documented default is not applied in 1.1.0.) Label -1 = unassigned.

Accuracy uses the pre-registered test-half pairs, weights and bootstrap of
v4_evaluate.py; stability uses benchmark_stability.compare.

Usage:
  python v4_mge_cluster.py
Output: output/backbone_v4/evaluation/{mge_cluster_test_metrics.tsv, stability_mge_cluster.tsv}
"""

import os

import numpy as np
import pandas as pd

from benchmark_stability import compare
from v4_evaluate import N_BOOT, SEED, TRUTHS, boot_f1, load_labels, metrics, same_group

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V4 = os.path.join(BASE_DIR, "output", "backbone_v4")
MGE = os.path.join(V4, "mge_cluster")


def labels(path):
    d = pd.read_csv(path, dtype=str)
    return {s: (None if c.strip() == "-1" else c.strip()) for s, c in zip(d.Sample_Name, d.Standard_Cluster)}


def main():
    pairs = pd.read_csv(os.path.join(V4, "truth_pairs.tsv"), sep="\t")
    test = pairs[pairs.half == "test"].reset_index(drop=True)
    a, b, w = test.plasmid_shorter, test.plasmid_longer, test.weight.values
    lab = load_labels(0.40)
    lab["mge-cluster"] = labels(os.path.join(MGE, "model_all", "mge-cluster_results.csv"))

    ids = sorted(set(a) | set(b))
    pos = {i: k for k, i in enumerate(ids)}
    rng = np.random.default_rng(SEED)
    counts = rng.multinomial(len(ids), np.full(len(ids), 1 / len(ids)), size=N_BOOT)
    mult = counts[:, [pos[x] for x in a]] * counts[:, [pos[x] for x in b]]

    rows, boots = [], {}
    for m in ("mge-cluster", "pLIN v4 L3", "pLIN v4 L6", "MOB-suite primary", "MOB-suite secondary"):
        pred = same_group(lab[m], a, b)
        for t in TRUTHS:
            y = test[t].values
            p, r, f = metrics(pred, y, w)
            fu = metrics(pred, y, np.ones_like(w))[2]
            boots[(m, t)] = boot_f1(pred, y, w, mult)
            rows.append({"method": m, "truth": t, "precision_w": p, "recall_w": r, "F1_w": f,
                         "F1_w_lo": np.percentile(boots[(m, t)], 2.5), "F1_w_hi": np.percentile(boots[(m, t)], 97.5),
                         "F1_unw": fu, "plasmids_assigned": np.mean([lab[m].get(i) is not None for i in ids])})
    for m1, t in (("pLIN v4 L3", "related_backbone"), ("pLIN v4 L6", "same_lineage")):
        d = boots[(m1, t)] - boots[("mge-cluster", t)]
        rows.append({"method": f"{m1} minus mge-cluster", "truth": t,
                     "F1_w": float(np.mean(d)), "F1_w_lo": np.percentile(d, 2.5), "F1_w_hi": np.percentile(d, 97.5)})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(V4, "evaluation", "mge_cluster_test_metrics.tsv"), sep="\t", index=False)

    A = [os.path.basename(l.strip())[:-len(".fasta")] for l in open(os.path.join(MGE, "snapshotA.txt")) if l.strip()]
    before = labels(os.path.join(MGE, "model_A", "mge-cluster_results.csv"))
    st = []
    for mode, f in (("existing", os.path.join(MGE, "existing_A_all", "mge-cluster_prediction.csv")),
                    ("rebuild", os.path.join(MGE, "model_all", "mge-cluster_results.csv"))):
        st.append({"tool": "mge-cluster", "mode": mode, **compare({i: before.get(i) for i in A}, labels(f))})
    st = pd.DataFrame(st)
    st.to_csv(os.path.join(V4, "evaluation", "stability_mge_cluster.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(res.round(3).to_string(index=False) + "\n\n" + st.to_string(index=False))


if __name__ == "__main__":
    main()
