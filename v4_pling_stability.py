#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
pling label stability on the pLIN v4 evaluation plasmids (secondary endpoint).

Snapshot A = random half of the 2,000 evaluation plasmids (seed 1). For the
plasmids in A, pling's labels from the run on A alone are compared with
  add      pling add (A + B on top of the A run; pling's update mode)
  rebuild  pling cluster on all 2,000 from scratch
using the same metrics as benchmark_stability.py.

Usage:
  python v4_pling_stability.py
Output: output/backbone_v4/evaluation/stability_pling.tsv
"""

import os

import pandas as pd

from benchmark_stability import compare, load_pling

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V4 = os.path.join(BASE_DIR, "output", "backbone_v4")
ST = os.path.join(V4, "stability")


def main():
    A = [os.path.basename(l.strip())[:-len(".fasta")] for l in open(os.path.join(ST, "snapshotA_seed1.txt")) if l.strip()]
    before = load_pling(os.path.join(ST, "pling_A"))
    rows = []
    for mode, d in (("add", os.path.join(ST, "pling_A_add")), ("rebuild", os.path.join(V4, "pling", "out"))):
        after = load_pling(d)
        for level in ("community", "subcommunity"):
            rows.append({"tool": "pling", "mode": mode, "level": level, "snapshot_frac": 0.5, "seed": 1,
                         **compare({i: before[level].get(i) for i in A}, after[level])})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(V4, "evaluation", "stability_pling.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
