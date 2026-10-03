#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Alignment truth set for the pLIN v4 evaluation (PREREGISTRATION.md, section 4).

Within each half of output/backbone_v4/eval_plasmids.tsv:
  1. every plasmid is sketched with sourmash (k = 21, scaled = 1000) and all
     pairs get a MinHash Jaccard similarity (a measure none of the compared
     tools uses for its final labels);
  2. pairs are stratified into Jaccard bins [0, 0.01), [0.01, 0.05),
     [0.05, 0.15), [0.15, 0.4), [0.4, 0.8), [0.8, 1] and up to 600 pairs per
     bin are drawn at random (seed 2026);
  3. each drawn pair is aligned exactly as in validate_alignment_backbone.py
     (blastn, one-to-one hit set; backbone masked where AMR/IS annotation
     exists).
Each pair carries weight = pairs in its bin / pairs drawn from its bin, so
weighted metrics estimate performance over all pairs in the half.

Truth (unchanged from earlier validation):
  same_lineage      AF_min >= 0.8 and ANI_aln >= 99
  related_backbone  bb_AF_shorter >= 0.5 (AF_shorter where no annotation)

Usage:
  python v4_truth_pairs.py --threads 12 [--per-bin 600] [--limit-plasmids N]
Output: output/backbone_v4/truth_pairs.tsv (+ sketches/, jaccard_<half>.csv)
"""

import argparse
import os
import random
import subprocess
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from Bio import SeqIO

from validate_alignment_backbone import align_pair, load_accessory_intervals

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "output", "backbone_v4")
SOURMASH = os.path.expanduser("~/miniconda3/envs/pling_env/bin/sourmash")
BIN_EDGES = [0, 0.01, 0.05, 0.15, 0.4, 0.8, 1.0000001]
SEED = 2026


def jaccard_matrix(paths, name, workdir):
    """All-pairs MinHash Jaccard for one half; returns (ids, matrix)."""
    lst = os.path.join(workdir, f"{name}_files.txt")
    sig = os.path.join(workdir, f"{name}.sig.zip")
    csv = os.path.join(workdir, f"jaccard_{name}.csv")
    open(lst, "w").write("\n".join(paths) + "\n")
    subprocess.run([SOURMASH, "sketch", "dna", "-p", "k=21,scaled=1000", "--from-file", lst, "-o", sig, "-f"],
                   check=True, capture_output=True)
    subprocess.run([SOURMASH, "compare", sig, "-k", "21", "--csv", csv], check=True, capture_output=True)
    m = pd.read_csv(csv)
    ids = [os.path.basename(c)[:-len(".fasta")] for c in m.columns]
    return ids, m.values


def sample_pairs(ids, mat, per_bin, rng):
    iu, ju = np.triu_indices(len(ids), k=1)
    jac = mat[iu, ju]
    b = np.digitize(jac, BIN_EDGES) - 1
    rows = []
    for k in range(len(BIN_EDGES) - 1):
        idx = np.flatnonzero(b == k)
        if not len(idx):
            continue
        pick = sorted(rng.sample(list(idx), min(per_bin, len(idx))))
        w = len(idx) / len(pick)
        for t in pick:
            rows.append({"p1": ids[iu[t]], "p2": ids[ju[t]], "jaccard": float(jac[t]),
                         "jaccard_bin": f"[{BIN_EDGES[k]}, {min(BIN_EDGES[k + 1], 1)})",
                         "bin_pairs": len(idx), "weight": w})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--per-bin", type=int, default=600)
    ap.add_argument("--limit-plasmids", type=int, default=0, help="testing only: first N plasmids per half")
    ap.add_argument("--out", default=os.path.join(OUT, "truth_pairs.tsv"))
    ap.add_argument("--eval", default=os.path.join(OUT, "eval_plasmids.tsv"),
                    help="plasmid table with columns plasmid_id, fasta, half")
    args = ap.parse_args()

    ev = pd.read_csv(args.eval, sep="\t")
    workdir = os.path.join(os.path.dirname(os.path.abspath(args.out)), "sketches")
    os.makedirs(workdir, exist_ok=True)
    rng = random.Random(SEED)
    pairs = []
    for half in [h for h in ("calibration", "test") if (ev.half == h).any()]:
        sub = ev[ev.half == half]
        if args.limit_plasmids:
            sub = sub.head(args.limit_plasmids)
        ids, mat = jaccard_matrix(sorted(sub.fasta), half + ("_test" if args.limit_plasmids else ""), workdir)
        assert set(ids) == set(sub.plasmid_id), "sourmash labels do not match plasmid IDs"
        for r in sample_pairs(ids, mat, args.per_bin, rng):
            pairs.append({"half": half, **r})
    pairs = pd.DataFrame(pairs)
    print(pairs.groupby(["half", "jaccard_bin"]).agg(drawn=("p1", "size"), in_bin=("bin_pairs", "first")).to_string())

    fasta = dict(zip(ev.plasmid_id, ev.fasta))
    length = {i: len(next(SeqIO.parse(fasta[i], "fasta")).seq) for i in set(pairs.p1) | set(pairs.p2)}
    intervals, _ = load_accessory_intervals()
    tasks, order = [], []
    for p1, p2 in zip(pairs.p1, pairs.p2):
        a, b = (p1, p2) if length[p1] <= length[p2] else (p2, p1)   # query = shorter plasmid
        order.append((a, b))
        tasks.append((fasta[a], fasta[b], length[a], length[b], intervals.get(a, []), intervals.get(b, []), False))
    with ProcessPoolExecutor(max_workers=args.threads) as ex:
        res = list(ex.map(align_pair, tasks, chunksize=4))

    out = []
    for (a, b), r, row in zip(order, res, pairs.itertuples(index=False)):
        out.append({"half": row.half, "plasmid_shorter": a, "plasmid_longer": b,
                    "jaccard": row.jaccard, "jaccard_bin": row.jaccard_bin, "weight": row.weight,
                    "len_shorter": length[a], "len_longer": length[b],
                    "AF_shorter": r["a"]["AF"], "AF_longer": r["b"]["AF"],
                    "AF_min": min(r["a"]["AF"], r["b"]["AF"]), "ANI_aln": r["ANI_aln"],
                    "bb_AF_shorter": r["a"]["bb_AF"], "annotated": bool(intervals.get(a) or intervals.get(b))})
    out = pd.DataFrame(out)
    out["same_lineage"] = (out.AF_min >= 0.8) & (out.ANI_aln >= 99)
    out["related_backbone"] = out.bb_AF_shorter.fillna(out.AF_shorter) >= 0.5
    out.to_csv(args.out, sep="\t", index=False)
    print(out.groupby("half")[["same_lineage", "related_backbone", "annotated"]].sum().to_string())
    print(f"wrote {len(out)} pairs to {args.out}")


if __name__ == "__main__":
    main()
