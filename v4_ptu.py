#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Agreement with published plasmid taxonomic units (PTUs; Redondo-Salvo et al.
2020, Nat Commun 11:3602, Supplementary Data 2).

  pre-registered  test-half evaluation plasmids that have a PTU (n = 36;
                  underpowered, reported as such)
  exploratory     every release plasmid with a PTU (n = 3,795); not
                  pre-registered. pling is not run at this scale.

For each method: adjusted Rand index against PTU, and pairwise precision /
recall / F1 for "same PTU" over all pairs of PTU-assigned plasmids
(unassigned labels count as singletons for ARI and never as "same group").

Usage:
  python v4_ptu.py
Output: output/backbone_v4/evaluation/ptu_agreement.tsv
"""

import glob
import os
from itertools import combinations

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V4 = os.path.join(BASE_DIR, "output", "backbone_v4")
PTU_XLSX = os.path.join(V4, "ptu", "sd_MOESM5.xlsx")
MOB_DIRS = [os.path.join(V4, "mob_typer"), os.path.join(V4, "ptu", "mob_typer"),
            os.path.join(BASE_DIR, "output", "comparator_benchmark", "mob_typer")]


def ptu_table():
    p = pd.read_excel(PTU_XLSX, header=1).iloc[:, :2]
    p.columns = ["acc", "PTU"]
    p = p.dropna(subset=["acc"])
    return p[p.PTU != "-"].set_index("acc").PTU


def release_ids_with_ptu():
    ptu = ptu_table()
    r = pd.read_csv(os.path.join(BASE_DIR, "output", "pLIN_reference_assignments.tsv"), sep="\t",
                    usecols=["plasmid_id"]).drop_duplicates()
    r["acc"] = r.plasmid_id.str.replace("^RefSeq_", "", regex=True)
    r = r[r.acc.isin(ptu.index)]
    return dict(zip(r.plasmid_id, r.acc.map(ptu)))


def mob_labels(ids):
    out = {}
    for i in ids:
        f = next((os.path.join(d, f"{i}.txt") for d in MOB_DIRS if os.path.exists(os.path.join(d, f"{i}.txt"))), None)
        if f is None:
            continue
        r = pd.read_csv(f, sep="\t", dtype=str).iloc[0]
        clean = lambda x: None if pd.isna(x) or x in ("-", "") else x
        out[i] = (clean(r.get("primary_cluster_id")), clean(r.get("secondary_cluster_id")))
    return out


def score(labels, truth, ids):
    a = [labels.get(i) if labels.get(i) is not None else f"_single_{i}" for i in ids]
    t = [truth[i] for i in ids]
    ari = adjusted_rand_score(t, a)
    tp = fp = fn = 0
    for x, y in combinations(range(len(ids)), 2):
        same_t = t[x] == t[y]
        same_a = labels.get(ids[x]) is not None and a[x] == a[y]
        tp += same_t and same_a
        fp += same_a and not same_t
        fn += same_t and not same_a
    p = tp / (tp + fp) if tp + fp else np.nan
    r = tp / (tp + fn) if tp + fn else np.nan
    return {"n": len(ids), "ARI": ari, "precision": p, "recall": r,
            "F1": 2 * p * r / (p + r) if tp else 0.0, "assigned_pct": 100 * np.mean([labels.get(i) is not None for i in ids])}


def methods(ids):
    lab = {}
    v3 = pd.read_csv(os.path.join(BASE_DIR, "output", "pLIN_reference_assignments.tsv"), sep="\t",
                     usecols=["plasmid_id", "pLIN"]).drop_duplicates("plasmid_id").set_index("plasmid_id").pLIN
    v4 = pd.read_csv(os.path.join(V4, "codes", "plin_v4_codes_L3_0.40.tsv"), sep="\t").set_index("plasmid_id").pLIN_v4
    for k in range(1, 7):
        lab[f"pLIN v3 L{k}"] = {i: ".".join(str(v3[i]).split(".")[:k]) for i in ids}
        lab[f"pLIN v4 L{k}"] = {i: ".".join(str(v4[i]).split(".")[:k]) for i in ids}
    mob = mob_labels(ids)
    lab["MOB-suite primary"] = {i: mob[i][0] for i in mob}
    lab["MOB-suite secondary"] = {i: mob[i][1] for i in mob}
    return lab, len(mob)


def main():
    truth = release_ids_with_ptu()
    ev = pd.read_csv(os.path.join(V4, "eval_plasmids.tsv"), sep="\t")
    sets = {"pre-registered: test half": sorted(set(ev[ev.half == "test"].plasmid_id) & set(truth)),
            "exploratory: all release plasmids with a PTU": sorted(truth)}
    rows = []
    for name, ids in sets.items():
        lab, n_mob = methods(ids)
        if n_mob < len(ids):
            print(f"{name}: MOB-suite output for {n_mob}/{len(ids)} plasmids; MOB-suite scored on those only")
        for m, l in lab.items():
            sub = [i for i in ids if not m.startswith("MOB") or i in l]
            rows.append({"set": name, "method": m, **score(l, truth, sub)})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(V4, "evaluation", "ptu_agreement.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
