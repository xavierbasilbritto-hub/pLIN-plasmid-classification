#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Discriminatory power of pLIN v4.1 levels and replicon typing (Simpson's index of diversity,
Hunter and Gaston 1988) over all unique plasmids of the release.

D = 1 - sum n_j (n_j - 1) / (N (N - 1)): the probability that two randomly chosen
plasmids get different types. Replicon type = the KNN replicon group in the release
table (plasmids without a replicon call count as one group, "Unknown").

Usage:
  python v41_diversity.py
Output: output/backbone_v41/diversity.json
"""

import json
import os

import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
REL = os.path.join(BASE_DIR, "output", "backbone_v41", "release", "plin_v41_codes.tsv.gz")


def simpson(labels):
    n = labels.value_counts().values.astype(float)
    N = n.sum()
    return 1 - (n * (n - 1)).sum() / (N * (N - 1))


def main():
    d = pd.read_csv(REL, sep="\t", dtype=str)
    out = {"plasmids": len(d), "replicon_typing": round(simpson(d.inc_type.fillna("Unknown")), 4),
           "replicon_types": int(d.inc_type.fillna("Unknown").nunique())}
    for k in range(1, 7):
        lab = d.pLIN_v41.str.split(".").str[:k].str.join(".")
        out[f"pLIN_L{k}"] = round(simpson(lab), 4)
        out[f"pLIN_L{k}_types"] = int(lab.nunique())
    json.dump(out, open(os.path.join(BASE_DIR, "output", "backbone_v41", "diversity.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
