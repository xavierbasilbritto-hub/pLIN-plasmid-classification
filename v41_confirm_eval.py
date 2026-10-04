#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
pLIN v4.1 confirmatory endpoints (PREREGISTRATION_v4.1.md, sections 6–8).

Test pairs: output/backbone_v41/confirm/truth_pairs.tsv (2,000 fresh plasmids).
Methods: pLIN v4.1 L1–L6 (codes_v41.tsv), pLIN v4 L1–L6 (pre-registered
baseline), MOB-suite primary/secondary, pling community/subcommunity,
mge-cluster (defaults). Weighted precision/recall/F1 with 95% bootstrap CIs
(1,000 resamples of plasmids, seed 2026) as in v4_evaluate.py.

Endpoints
  E1 related backbone: F1 v4.1 L1 − F1 MOB-suite primary   (R5: lower bound > −0.02)
  E2 same lineage:     F1 v4.1 L5 − F1 MOB-suite secondary (non-inferiority reported)

Usage:
  python v41_confirm_eval.py
Output: output/backbone_v41/confirm/{test_metrics.tsv, endpoints.json}
"""

import glob
import json
import os

import numpy as np
import pandas as pd

from v4_evaluate import N_BOOT, SEED, TRUTHS, boot_f1, metrics, same_group

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V41 = os.path.join(BASE_DIR, "output", "backbone_v41")
C = os.path.join(V41, "confirm")
MARGIN = 0.02


def prefix(code, k):
    return ".".join(str(code).split(".")[:k])


def labels(ids):
    lab = {}
    v41 = pd.read_csv(os.path.join(V41, "codes_v41.tsv"), sep="\t").set_index("plasmid_id").pLIN_v41
    v4 = pd.read_csv(os.path.join(BASE_DIR, "output", "backbone_v4", "codes", "plin_v4_codes_L3_0.40.tsv"),
                     sep="\t").set_index("plasmid_id").pLIN_v4
    for k in range(1, 7):
        lab[f"pLIN v4.1 L{k}"] = {i: prefix(v41[i], k) for i in ids}
        lab[f"pLIN v4 L{k}"] = {i: prefix(v4[i], k) for i in ids}
    clean = lambda x: None if pd.isna(x) or x in ("-", "") else x
    mob = {}
    for i in ids:
        r = pd.read_csv(os.path.join(C, "mob_typer", f"{i}.txt"), sep="\t", dtype=str).iloc[0]
        mob[i] = (clean(r.get("primary_cluster_id")), clean(r.get("secondary_cluster_id")))
    lab["MOB-suite primary"] = {i: mob[i][0] for i in ids}
    lab["MOB-suite secondary"] = {i: mob[i][1] for i in ids}
    out = os.path.join(C, "pling", "out")
    typing = pd.read_csv(glob.glob(os.path.join(out, "dcj_thresh_*_graph", "objects", "typing.tsv"))[0], sep="\t", dtype=str)
    comm = pd.read_csv(os.path.join(out, "containment", "containment_communities", "objects", "communities.tsv"),
                       sep="\t", dtype=str)
    sub, com = dict(typing.iloc[:, :2].values), dict(comm.iloc[:, :2].values)
    lab["pling subcommunity"] = {i: sub.get(i) for i in ids}
    lab["pling community"] = {i: com.get(i) for i in ids}
    m = pd.read_csv(os.path.join(C, "mge", "model_all", "mge-cluster_results.csv"), dtype=str)
    mg = {s: (None if c.strip() == "-1" else c.strip()) for s, c in zip(m.Sample_Name, m.Standard_Cluster)}
    lab["mge-cluster"] = {i: mg.get(i) for i in ids}
    return lab


def main():
    pairs = pd.read_csv(os.path.join(C, "truth_pairs.tsv"), sep="\t").reset_index(drop=True)
    a, b, w = pairs.plasmid_shorter, pairs.plasmid_longer, pairs.weight.values
    ids = sorted(set(a) | set(b))
    lab = labels(sorted(pd.read_csv(os.path.join(C, "test_plasmids.tsv"), sep="\t").plasmid_id))
    pos = {i: k for k, i in enumerate(ids)}
    rng = np.random.default_rng(SEED)
    counts = rng.multinomial(len(ids), np.full(len(ids), 1 / len(ids)), size=N_BOOT)
    mult = counts[:, [pos[x] for x in a]] * counts[:, [pos[x] for x in b]]

    rows, boots = [], {}
    for name, l in lab.items():
        pred = same_group(l, a, b)
        for t in TRUTHS:
            y = pairs[t].values
            p, r, f = metrics(pred, y, w)
            pu, ru, fu = metrics(pred, y, np.ones_like(w))
            boots[(name, t)] = boot_f1(pred, y, w, mult)
            rows.append({"method": name, "truth": t, "precision_w": p, "recall_w": r, "F1_w": f,
                         "F1_w_lo": np.percentile(boots[(name, t)], 2.5), "F1_w_hi": np.percentile(boots[(name, t)], 97.5),
                         "precision_unw": pu, "recall_unw": ru, "F1_unw": fu,
                         "plasmids_assigned": float(np.mean([l.get(i) is not None for i in ids]))})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(C, "test_metrics.tsv"), sep="\t", index=False)

    def endpoint(m1, m2, t):
        get = lambda m: float(res[(res.method == m) & (res.truth == t)].F1_w.iloc[0])
        d = boots[(m1, t)] - boots[(m2, t)]
        return {"F1": get(m1), "F1_comparator": get(m2), "difference": get(m1) - get(m2),
                "lo": float(np.percentile(d, 2.5)), "hi": float(np.percentile(d, 97.5)),
                "noninferior": bool(np.percentile(d, 2.5) > -MARGIN)}
    ep = {"test_pairs": int(len(pairs)), "same_lineage_pairs": int(pairs.same_lineage.sum()),
          "related_backbone_pairs": int(pairs.related_backbone.sum()),
          "E1_backbone_v41L1_vs_MOBprimary": endpoint("pLIN v4.1 L1", "MOB-suite primary", "related_backbone"),
          "E2_lineage_v41L5_vs_MOBsecondary": endpoint("pLIN v4.1 L5", "MOB-suite secondary", "same_lineage")}
    ep["R5_pass"] = ep["E1_backbone_v41L1_vs_MOBprimary"]["noninferior"]
    json.dump(ep, open(os.path.join(C, "endpoints.json"), "w"), indent=2)
    pd.set_option("display.width", 250)
    print(res.round(3).to_string(index=False))
    print(json.dumps(ep, indent=2))


if __name__ == "__main__":
    main()
