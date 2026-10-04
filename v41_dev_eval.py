#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Evaluate one candidate v4.1 design on DEVELOPMENT data (all previously seen).
Engine: "founder" (plin_v41.V41Tree) or "nn" (plin_v41_nn.NNIndex, LIN-code style).

Builds codes for every unique release plasmid (accession order), then reports
  pairs      weighted F1 per level vs alignment truth (5,544 v4 evaluation pairs)
  requery    share of 2,000 release plasmids that re-query to their own code
  mutants    share of mutated plasmids keeping the original's code, per level
  vim        Swiss VIM-1: for the 8 blaVIM-1 plasmids, does the deepest shared
             level match the alignment (same lineage vs related only)?
  outbreak   same-study / different-study pairs grouped, per level
  clusters   number of clusters per level

Usage:
  python v41_dev_eval.py '<json config>' NAME
  config: {"levels": [["prot", 0.5], ..., ["kmin", 0.9]], "rule": "best"|"earliest", "engine": "founder"|"nn"}
Output: output/backbone_v41/dev/eval_NAME.json
"""

import json
import os
import random
import sys
import time
from itertools import combinations

import numpy as np
import pandas as pd

import plin_v4 as v
from plin_kmers import containment, min_containment
from plin_v41 import V41Tree
from plin_v41_nn import MIN_SHARED, NNIndex

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DEV = os.path.join(BASE_DIR, "output", "backbone_v41", "dev")
V4 = os.path.join(BASE_DIR, "output", "backbone_v4")


def acc(i):
    return i.replace("RefSeq_", "", 1)


def load_db():
    z = np.load(os.path.join(BASE_DIR, "output", "backbone_v41", "sketches.npz"))
    ids, off, h, sc = z["ids"].tolist(), z["offsets"], z["hashes"], z["scaled"]
    sk = {i: (h[off[k]:off[k + 1]], int(sc[k])) for k, i in enumerate(ids)}
    fs = v.load_family_sets()
    fams = {i: fs[i] for i in ids}
    return ids, fams, sk


def load_queries():
    z = np.load(os.path.join(DEV, "query_inputs.npz"))
    q = {}
    for k, (g, r) in enumerate(zip(z["group"], z["rid"])):
        q[str(r)] = {"group": str(g), "fams": z["fam"][z["fam_off"][k]:z["fam_off"][k + 1]].astype(np.int64),
                     "sk": z["sk"][z["sk_off"][k]:z["sk_off"][k + 1]], "scale": int(z["scale"][k])}
    return q


def shared_levels(a, b):
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


def f1(pred, truth, w):
    tp, pp, pos = (w * (pred & truth)).sum(), (w * pred).sum(), (w * truth).sum()
    p, r = (tp / pp if pp else 0.0), (tp / pos if pos else 0.0)
    return round(2 * p * r / (p + r), 4) if p + r else 0.0, round(p, 4), round(r, 4)


def main():
    cfg, name = json.loads(sys.argv[1]), sys.argv[2]
    ids, fams, sk = load_db()
    t0 = time.perf_counter()
    levels = [tuple(x) for x in cfg["levels"]]
    if cfg.get("engine") == "hybrid":
        from plin_v41_nn import HybridIndex
        F, S, C = [fams[i] for i in ids], [sk[i][0] for i in ids], [sk[i][1] for i in ids]
        hx = HybridIndex(levels, cfg["lin"], F, S, C)
        hx.build(F, S, C, progress=20000)
        codes = {i: tuple(int(x) for x in hx.codes[k]) for k, i in enumerate(ids)}
        hsession = {"sk": [], "scale": [], "codes": []}

        def assign(f, s, c, add=True):
            extra = None
            if add and hsession["codes"]:
                extra = {"codes": np.array(hsession["codes"]),
                         "kmin": np.array([min_containment(s, s2, c, c2) for s2, c2 in zip(hsession["sk"], hsession["scale"])])}
            code, _ = hx.assign(f, s, c, len(F), extra=extra, add=add)
            if add:
                hsession["sk"].append(s); hsession["scale"].append(c); hsession["codes"].append(code)
            return tuple(code)
        levels = levels + [("kmin", t) for t in cfg["lin"]]
    elif cfg.get("engine") == "nn":
        F, S, C = [fams[i] for i in ids], [sk[i][0] for i in ids], [sk[i][1] for i in ids]
        ix = NNIndex(levels, F, S, C)
        ix.build(F, S, C, progress=20000)
        codes = {i: tuple(int(x) for x in ix.codes[k]) for k, i in enumerate(ids)}
        pos = {i: k for k, i in enumerate(ids)}
        session = {"fams": [], "sk": [], "scale": [], "codes": []}

        def pair_sims(f, s, c):
            if not session["codes"]:
                return None
            out = {"prot": [], "kcont": [], "kmin": []}
            for f2, s2, c2 in zip(session["fams"], session["sk"], session["scale"]):
                inter = len(np.intersect1d(f, f2))
                small = min(len(f), len(f2))
                out["prot"].append(inter / small if small and inter >= min(MIN_SHARED, small) else 0.0)
                out["kcont"].append(containment(s, s2, c, c2) if len(s) <= len(s2) else containment(s2, s, c2, c))
                out["kmin"].append(min_containment(s, s2, c, c2))
            return {"codes": np.array(session["codes"]), "sims": {k: np.array(x) for k, x in out.items()}}

        def assign(f, s, c, add=True):
            sims = ix.similarities(f, s, c, ix.n)
            extra = pair_sims(f, s, c) if add else None
            code, _ = ix.code_for(sims, ix.n, noprot=len(f) == 0, extra=extra,
                                  next_id=ix.next_id if add else ix.next_id.copy())
            if add:
                session["fams"].append(f); session["sk"].append(s); session["scale"].append(c)
                session["codes"].append(code)
            return tuple(code)
    else:
        tree = V41Tree(levels, rule=cfg.get("rule", "best"))
        codes = {i: tree.assign(fams[i], *sk[i]) for i in ids}

        def assign(f, s, c, add=True):
            return tree.assign(f, s, c, add=add)
    build_s = time.perf_counter() - t0
    nlev = len(levels)
    by_acc = {acc(i): codes[i] for i in ids}
    res = {"config": cfg, "build_seconds": round(build_s, 1),
           "clusters": {f"L{k + 1}": len({c[:k + 1] for c in codes.values()}) for k in range(nlev)}}

    pairs = pd.read_csv(os.path.join(V4, "truth_pairs.tsv"), sep="\t")
    ca = [by_acc[acc(x)] for x in pairs.plasmid_shorter]
    cb = [by_acc[acc(x)] for x in pairs.plasmid_longer]
    sh = np.array([shared_levels(a, b) for a, b in zip(ca, cb)])
    res["pairs"] = {}
    for half in ("calibration", "test", "all"):
        m = (pairs.half == half).values if half != "all" else np.ones(len(pairs), bool)
        w = pairs.weight.values[m]
        res["pairs"][half] = {t: {f"L{k + 1}": f1(sh[m] >= k + 1, pairs[t].values[m], w) for k in range(nlev)}
                              for t in ("related_backbone", "same_lineage")}

    sample = random.Random(1).sample(ids, 2000)
    res["requery_identical"] = float(np.mean([assign(fams[i], *sk[i], add=False) == codes[i] for i in sample]))

    q = load_queries()
    res["mutants"] = {}
    for rate in ("0.0001", "0.001", "0.01"):
        keep = []
        for r, d in q.items():
            if d["group"] == f"mut_{rate}":
                o = codes[r.split("__")[0]]
                c = assign(d["fams"], d["sk"], d["scale"], add=False)
                keep.append([c[:k + 1] == o[:k + 1] for k in range(nlev)])
        res["mutants"][rate] = [round(x, 3) for x in np.mean(keep, axis=0)]

    qc = {}
    for grp in ("outbreak", "vim"):                          # session: queries extend the tree
        for r, d in sorted(q.items()):
            if d["group"] == grp:
                qc[r] = assign(d["fams"], d["sk"], d["scale"], add=True)

    vp = pd.read_csv(os.path.join(V4, "swiss_vim1", "vim_pairs.tsv"), sep="\t")
    vrows = [{"a": a, "b": b, "same_lineage": bool(s), "shared": shared_levels(qc[a], qc[b])}
             for a, b, s in zip(vp.contig_shorter, vp.contig_longer, vp.same_lineage_by_alignment)]
    vd = pd.DataFrame(vrows)
    res["vim"] = {"same_lineage_pairs": int(vd.same_lineage.sum()),
                  "same_lineage_pairs_sharing_all_levels": int((vd[vd.same_lineage].shared == nlev).sum()),
                  "other_pairs": int((~vd.same_lineage).sum()),
                  "other_pairs_split_at_last_level": int((vd[~vd.same_lineage].shared < nlev).sum()),
                  "all_vim_pairs_min_shared_level": int(vd.shared.min()),
                  "shared_levels_by_pair": vd.shared.tolist()}

    ob = pd.read_csv(os.path.join(BASE_DIR, "output", "outbreak_validation_founder_results.tsv"), sep="\t")
    st = dict(zip(ob.accession, ob.study))
    res["outbreak"] = {}
    for k in range(nlev):
        same, diff = [], []
        for a, b in combinations(ob.accession, 2):
            tog = qc[a][:k + 1] == qc[b][:k + 1]
            (same if st[a] == st[b] else diff).append(tog)
        res["outbreak"][f"L{k + 1}"] = [round(100 * np.mean(same), 1), round(100 * np.mean(diff), 2)]

    json.dump(res, open(os.path.join(DEV, f"eval_{name}.json"), "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("build_seconds", "clusters", "requery_identical", "mutants", "vim", "outbreak")}, indent=None))
    for t in ("related_backbone", "same_lineage"):
        print(t, "F1 (all dev pairs):", {k: v_[0] for k, v_ in res["pairs"]["all"][t].items()})


if __name__ == "__main__":
    main()
