#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
pLIN v4.1 release criteria R1–R3 on the confirmatory data (PREREGISTRATION_v4.1.md, section 8).

R1 stability   the 2,000 test plasmids: for snapshots of 25/50/75% (seeds 1–5) the
               snapshot is coded first, the rest is added after it, and every snapshot
               plasmid must re-query to the code it had before the growth (all levels).
               A rebuild in accession order is reported for reference.
R2 reproduction 500 fresh release plasmids typed from sequence (pyrodigal + catalogue
               families + k-mer sketch; no whole-plasmid shortcut) must get their own code.
R3 robustness  100 fresh plasmids with random substitutions at 0.01%, 0.1%, 1%
               (seed 4343): share keeping the original code at L1..L6; criterion
               >= 95% for L1–L5 at 0.01% and 0.1%.

Usage:
  python v41_confirm_checks.py
Output: output/backbone_v41/confirm/release_checks.json
"""

import json
import os
import pickle
import random

import numpy as np
import pandas as pd
from Bio import SeqIO

import plin_v4 as v
from plin_kmers import adaptive_sketch
from plin_v41_nn import HybridIndex
from plin_v4_typer import PlinV4Release
from v41_build import FOUNDER_LEVELS, LINEAGE

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V41 = os.path.join(BASE_DIR, "output", "backbone_v41")
C = os.path.join(V41, "confirm")


def mutate(s, rate, r):
    s = list(s)
    for k in r.sample(range(len(s)), max(1, int(len(s) * rate))):
        s[k] = r.choice([b for b in "ACGT" if b != s[k]])
    return "".join(s)


def main():
    with open(os.path.join(V41, "hybrid_index.pkl"), "rb") as fh:
        d = pickle.load(fh)
    ids, hx = d["ids"], d["index"]
    pos = {i: k for k, i in enumerate(ids)}
    n = len(ids)
    res = {}

    # R1: stability on the 2,000 test plasmids
    z = np.load(os.path.join(V41, "sketches.npz"))
    off, h, sc = z["offsets"], z["hashes"], z["scaled"]
    fs = v.load_family_sets()
    test = pd.read_csv(os.path.join(C, "test_plasmids.tsv"), sep="\t").plasmid_id.tolist()
    get = lambda i: (fs[i], h[off[pos[i]]:off[pos[i] + 1]], int(sc[pos[i]]))
    rows = []
    for frac in (0.25, 0.5, 0.75):
        for seed in (1, 2, 3, 4, 5):
            A = sorted(random.Random(seed).sample(test, int(len(test) * frac)))
            B = sorted(set(test) - set(A))
            order = A + B
            F, S, Cc = zip(*[get(i) for i in order])
            g = HybridIndex(FOUNDER_LEVELS, LINEAGE, list(F), list(S), list(Cc))
            g.build(list(F), list(S), list(Cc))
            before = g.codes[:len(A)].copy()                 # codes of A, coded before B existed
            after = np.array([g.assign(F[k], S[k], Cc[k], len(order), add=False)[0] for k in range(len(A))])
            rows.append({"frac": frac, "seed": seed,
                         "retained_by_level": [float((before[:, :l + 1] == after[:, :l + 1]).all(axis=1).mean())
                                               for l in range(6)]})
    res["R1_stability"] = {"runs": rows,
                           "min_retained_all_levels": min(r["retained_by_level"][5] for r in rows),
                           "pass": all(r["retained_by_level"][5] == 1.0 for r in rows)}
    print("R1", res["R1_stability"]["min_retained_all_levels"], flush=True)

    # R2 / R3: sequence path against the full release index
    rel = PlinV4Release(os.path.join(BASE_DIR, "output", "backbone_v4", "release"))
    codes = {i: tuple(int(x) for x in hx.codes[pos[i]]) for i in ids}

    def type_seqs(records, work):
        os.makedirs(work, exist_ok=True)
        paths = {}
        for k, (rid, seq) in enumerate(records):
            p = os.path.join(work, f"q{k}.fa")
            open(p, "w").write(f">{rid}\n{seq}\n")
            paths[rid] = p
        fams = rel._families(paths, work, threads=8)
        out = {}
        for rid, seq in records:
            sk, s = adaptive_sketch(seq)
            out[rid] = tuple(hx.assign(fams[rid], sk, s, n, add=False)[0])
        return out

    req = pd.read_csv(os.path.join(C, "requery_plasmids.tsv"), sep="\t")
    recs = [(i, "".join(str(r.seq) for r in SeqIO.parse(f, "fasta")).upper()) for i, f in zip(req.plasmid_id, req.fasta)]
    got = type_seqs(recs, os.path.join(C, "work_requery"))
    same = [got[i] == codes[i] for i, _ in recs]
    res["R2_reproduction"] = {"n": len(same), "identical": float(np.mean(same)), "pass": all(same),
                              "mismatches": [i for (i, _), s in zip(recs, same) if not s][:20]}
    print("R2", res["R2_reproduction"]["identical"], flush=True)

    mut = pd.read_csv(os.path.join(C, "mutation_plasmids.tsv"), sep="\t")
    orig = {i: "".join(str(r.seq) for r in SeqIO.parse(f, "fasta")).upper() for i, f in zip(mut.plasmid_id, mut.fasta)}
    res["R3_robustness"] = {}
    for rate in (0.0001, 0.001, 0.01):
        r = random.Random(4343)
        recs = [(f"{i}__{rate}", mutate(orig[i], rate, r)) for i in sorted(orig)]
        got = type_seqs(recs, os.path.join(C, f"work_mut_{rate}"))
        keep = [[got[q][:l + 1] == codes[q.split("__")[0]][:l + 1] for l in range(6)] for q, _ in recs]
        res["R3_robustness"][str(rate)] = [round(float(x), 3) for x in np.mean(keep, axis=0)]
        print("R3", rate, res["R3_robustness"][str(rate)], flush=True)
    r3 = res["R3_robustness"]
    res["R3_pass"] = all(min(r3[k][:5]) >= 0.95 for k in ("0.0001", "0.001"))
    json.dump(res, open(os.path.join(C, "release_checks.json"), "w"), indent=1)
    print(json.dumps({k: (v_ if not isinstance(v_, dict) or "runs" not in v_ else {kk: vv for kk, vv in v_.items() if kk != "runs"})
                      for k, v_ in res.items()}, indent=1))


if __name__ == "__main__":
    main()
