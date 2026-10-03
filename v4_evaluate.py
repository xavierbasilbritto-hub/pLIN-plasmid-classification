#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Pre-registered evaluation of pLIN v4 (output/backbone_v4/PREREGISTRATION.md).

  1. Calibration half only: choose L3 from the pre-registered grid by weighted
     F1 against "related backbone" (ties -> lowest threshold).
  2. Test half: pairwise precision / recall / F1 for every method and both
     truths, weighted by inverse sampling probability, with 95% CIs from
     1,000 bootstrap resamples of test plasmids (a pair's weight is multiplied
     by how often each of its plasmids is drawn).
  3. Label stability of v4 on the 2,000 evaluation plasmids (design of
     benchmark_stability.py: snapshots 25/50/75%, seeds 1-5; incremental and
     rebuilt).
  4. Gate 1, applied mechanically: lower 95% bound of
     F1(v4 L3) - F1(MOB-suite primary) on related backbone > -0.02, and v4
     incremental stability 100% at every level.

Labels: pLIN v3 from output/pLIN_reference_assignments.tsv (release
db-2026.10.02); pLIN v4 from output/backbone_v4/codes/; MOB-suite mob_typer
(output/backbone_v4/mob_typer); pling (output/backbone_v4/pling/out; hub
plasmids, which pling leaves untyped, count as unassigned).

Usage:
  python v4_evaluate.py
Output: output/backbone_v4/evaluation/{calibration.tsv, test_metrics.tsv,
        primary_endpoints.json, stability_v4.tsv}
"""

import glob
import json
import os
import random
import warnings

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from plin_v4 import L3_GRID, V4Tree, load_family_sets, load_vectors

warnings.filterwarnings("ignore", message="The number of unique classes")
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V4 = os.path.join(BASE_DIR, "output", "backbone_v4")
EV = os.path.join(V4, "evaluation")
N_BOOT, MARGIN, SEED = 1000, 0.02, 2026
TRUTHS = ["related_backbone", "same_lineage"]


def prefix(code, k):
    return ".".join(str(code).split(".")[:k])


def load_labels(l3):
    ev = pd.read_csv(os.path.join(V4, "eval_plasmids.tsv"), sep="\t")
    ids = set(ev.plasmid_id)
    lab = {}
    v3 = pd.read_csv(os.path.join(BASE_DIR, "output", "pLIN_reference_assignments.tsv"), sep="\t",
                     usecols=["plasmid_id", "pLIN"]).drop_duplicates("plasmid_id").set_index("plasmid_id").pLIN
    for k in (3, 4, 5, 6):
        lab[f"pLIN v3 L{k}"] = {i: prefix(v3[i], k) for i in ids}
    v4 = pd.read_csv(os.path.join(V4, "codes", f"plin_v4_codes_L3_{l3:.2f}.tsv"), sep="\t").set_index("plasmid_id").pLIN_v4
    for k in range(1, 7):
        lab[f"pLIN v4 L{k}"] = {i: prefix(v4[i], k) for i in ids}
    mob = {}
    for i in ids:
        r = pd.read_csv(os.path.join(V4, "mob_typer", f"{i}.txt"), sep="\t", dtype=str).iloc[0]
        mob[i] = (r.get("primary_cluster_id"), r.get("secondary_cluster_id"))
    clean = lambda x: None if pd.isna(x) or x in ("-", "") else x
    lab["MOB-suite primary"] = {i: clean(mob[i][0]) for i in ids}
    lab["MOB-suite secondary"] = {i: clean(mob[i][1]) for i in ids}
    out = os.path.join(V4, "pling", "out")
    typing = pd.read_csv(glob.glob(os.path.join(out, "dcj_thresh_*_graph", "objects", "typing.tsv"))[0], sep="\t", dtype=str)
    comm = pd.read_csv(os.path.join(out, "containment", "containment_communities", "objects", "communities.tsv"),
                       sep="\t", dtype=str)
    sub, com = dict(typing.iloc[:, :2].values), dict(comm.iloc[:, :2].values)
    lab["pling subcommunity"] = {i: sub.get(i) for i in ids}
    lab["pling community"] = {i: com.get(i) for i in ids}
    return lab


def same_group(lab, a, b):
    return np.array([lab[x] is not None and lab[x] == lab[y] for x, y in zip(a, b)])


def metrics(pred, truth, w):
    tp, pp, pos = (w * (pred & truth)).sum(), (w * pred).sum(), (w * truth).sum()
    p = tp / pp if pp else np.nan
    r = tp / pos if pos else np.nan
    f = 2 * p * r / (p + r) if pp and pos and (p + r) else 0.0
    return p, r, f


def boot_f1(pred, truth, w, mult):
    """F1 per bootstrap replicate; mult is (R, n_pairs) pair multiplicities."""
    W = mult * w
    tp, pp, pos = (W * (pred & truth)).sum(1), (W * pred).sum(1), (W * truth).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        p, r = tp / pp, tp / pos
        f = np.where((p + r) > 0, 2 * p * r / (p + r), 0.0)
    return np.nan_to_num(f)


def calibrate(pairs):
    cal = pairs[pairs.half == "calibration"]
    rows = []
    for t in L3_GRID:
        lab = load_labels(t)["pLIN v4 L3"]
        p, r, f = metrics(same_group(lab, cal.plasmid_shorter, cal.plasmid_longer),
                          cal.related_backbone.values, cal.weight.values)
        rows.append({"L3": t, "precision": p, "recall": r, "F1": f})
    df = pd.DataFrame(rows)
    best = float(df.loc[df.F1.idxmax(), "L3"])            # idxmax -> first (lowest) on ties
    return df, best


def evaluate_test(pairs, l3):
    test = pairs[pairs.half == "test"].reset_index(drop=True)
    lab = load_labels(l3)
    a, b, w = test.plasmid_shorter, test.plasmid_longer, test.weight.values
    preds = {m: same_group(l, a, b) for m, l in lab.items()}
    preds["pLIN v4 L6 AND MOB-suite secondary"] = preds["pLIN v4 L6"] & preds["MOB-suite secondary"]
    preds["pLIN v3 L6 AND MOB-suite secondary"] = preds["pLIN v3 L6"] & preds["MOB-suite secondary"]

    test_ids = sorted(set(a) | set(b))
    pos = {i: k for k, i in enumerate(test_ids)}
    ia, ib = np.array([pos[x] for x in a]), np.array([pos[x] for x in b])
    rng = np.random.default_rng(SEED)
    counts = rng.multinomial(len(test_ids), np.full(len(test_ids), 1 / len(test_ids)), size=N_BOOT)
    mult = counts[:, ia] * counts[:, ib]

    rows, boots = [], {}
    for m, pred in preds.items():
        cov = np.mean([lab[m][i] is not None for i in test_ids]) if m in lab else np.nan
        for t in TRUTHS:
            truth = test[t].values
            p, r, f = metrics(pred, truth, w)
            pu, ru, fu = metrics(pred, truth, np.ones_like(w))
            bf = boot_f1(pred, truth, w, mult)
            boots[(m, t)] = bf
            rows.append({"method": m, "truth": t, "precision_w": p, "recall_w": r, "F1_w": f,
                         "F1_w_lo": np.percentile(bf, 2.5), "F1_w_hi": np.percentile(bf, 97.5),
                         "precision_unw": pu, "recall_unw": ru, "F1_unw": fu,
                         "plasmids_assigned": cov})
    return pd.DataFrame(rows), boots, len(test)


def diff_ci(boots, m1, m2, truth, point):
    d = boots[(m1, truth)] - boots[(m2, truth)]
    return {"difference": point, "lo": float(np.percentile(d, 2.5)), "hi": float(np.percentile(d, 97.5))}


def stability_v4(l3):
    """Incremental vs rebuilt v4 codes on the 2,000 evaluation plasmids."""
    ev = pd.read_csv(os.path.join(V4, "eval_plasmids.tsv"), sep="\t")
    ids_all = sorted(ev.plasmid_id)
    fams_all, vecs_all = load_family_sets(), load_vectors()
    fams = {i: fams_all[i] for i in ids_all}
    vecs = {i: vecs_all[i] for i in ids_all}

    def build(ids):
        tree = V4Tree(l3)
        return {i: tree.assign(fams[i], vecs[i]) for i in ids}, tree

    rebuilt, _ = build(ids_all)
    rows = []
    for frac in (0.25, 0.5, 0.75):
        for seed in (1, 2, 3, 4, 5):
            A = sorted(random.Random(seed).sample(ids_all, int(len(ids_all) * frac)))
            codes_a, tree = build(A)
            for i in sorted(set(ids_all) - set(A)):
                tree.assign(fams[i], vecs[i])
            requery = {i: tree.assign(fams[i], vecs[i], add=False) for i in A}
            for mode, after in (("incremental", requery), ("rebuild", rebuilt)):
                for k in range(1, 7):
                    before_l = [codes_a[i][:k] for i in A]
                    after_l = [after[i][:k] for i in A]
                    rows.append({"mode": mode, "level": f"L{k}", "snapshot_frac": frac, "seed": seed,
                                 "label_retained_pct": 100 * np.mean([x == y for x, y in zip(before_l, after_l)]),
                                 "ARI": adjusted_rand_score([str(x) for x in before_l], [str(x) for x in after_l])})
    return pd.DataFrame(rows)


def main():
    os.makedirs(EV, exist_ok=True)
    pairs = pd.read_csv(os.path.join(V4, "truth_pairs.tsv"), sep="\t")

    cal, l3 = calibrate(pairs)
    cal.to_csv(os.path.join(EV, "calibration.tsv"), sep="\t", index=False)
    print("Calibration half (choose L3):\n" + cal.round(4).to_string(index=False) + f"\n-> L3 = {l3:.2f}\n")

    res, boots, n_pairs = evaluate_test(pairs, l3)
    res.to_csv(os.path.join(EV, "test_metrics.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 250)
    print("Test half:\n" + res.round(4).to_string(index=False) + "\n")

    get = lambda m, t: float(res[(res.method == m) & (res.truth == t)].F1_w.iloc[0])
    e1 = diff_ci(boots, "pLIN v4 L3", "MOB-suite primary", "related_backbone",
                 get("pLIN v4 L3", "related_backbone") - get("MOB-suite primary", "related_backbone"))
    e2 = diff_ci(boots, "pLIN v4 L6", "MOB-suite secondary", "same_lineage",
                 get("pLIN v4 L6", "same_lineage") - get("MOB-suite secondary", "same_lineage"))

    st = stability_v4(l3)
    st.to_csv(os.path.join(EV, "stability_v4.tsv"), sep="\t", index=False)
    inc = st[st["mode"] == "incremental"]
    stable = bool((inc.label_retained_pct == 100).all())
    print("v4 stability (min over seeds):\n" +
          st.groupby(["mode", "level"])[["label_retained_pct", "ARI"]].min().round(3).to_string() + "\n")

    gate1 = e1["lo"] > -MARGIN and stable
    summary = {"L3_selected": l3, "test_pairs": n_pairs,
               "endpoint1_related_backbone_F1_v4L3_minus_MOBprimary": e1,
               "endpoint2_same_lineage_F1_v4L6_minus_MOBsecondary": e2,
               "v4_incremental_stability_100pct": stable, "margin": MARGIN,
               "gate1_pass": gate1, "endpoint2_noninferior": e2["lo"] > -MARGIN}
    json.dump(summary, open(os.path.join(EV, "primary_endpoints.json"), "w"), indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
