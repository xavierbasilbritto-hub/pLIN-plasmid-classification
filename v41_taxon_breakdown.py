#!/usr/bin/env python3
# Copyright (C) 2025-2026 Basil Britto Xavier. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Lineage accuracy by organism group (secondary analysis, not pre-registered).

The confirmatory test set is drawn from the whole release, so it is not limited to
Enterobacterales. This splits the registered same-lineage comparison by the replicon group of
the two plasmids of each pair, to show whether the method works outside Gram-negative plasmids.

Usage:
  python v41_taxon_breakdown.py
Output: output/backbone_v41/confirm/taxon_breakdown.json
"""

import json
import os

import numpy as np
import pandas as pd

from v41_confirm_eval import C, labels
from v4_evaluate import metrics, same_group

GRAM_POSITIVE = {"repSA_large", "repSA_small", "repEF_conj", "repEF_res"}
NONFERMENTER = {"repAci1", "repAci_large", "repPae_large", "repPae_small"}
BASE = os.path.dirname(os.path.abspath(__file__))


def main():
    pairs = pd.read_csv(os.path.join(C, "truth_pairs.tsv"), sep="\t")
    rel = pd.read_csv(os.path.join(BASE, "output", "backbone_v41", "release", "plin_v41_codes.tsv.gz"),
                      sep="\t", usecols=["accession", "inc_type"])
    inc = dict(zip(rel.accession, rel.inc_type))
    lab = labels(sorted(pd.read_csv(os.path.join(C, "test_plasmids.tsv"), sep="\t").plasmid_id))
    a, b, w = pairs.plasmid_shorter, pairs.plasmid_longer, pairs.weight.values

    def group(x):
        t = inc.get(x)
        return "gram_positive" if t in GRAM_POSITIVE else ("nonfermenter" if t in NONFERMENTER else "enterobacterales")

    out = {}
    for name in ("gram_positive", "nonfermenter", "enterobacterales"):
        m = np.array([group(x) == name and group(y) == name for x, y in zip(a, b)])
        if m.sum() < 30:
            continue
        row = {"pairs": int(m.sum()), "same_lineage_pairs": int(pairs.same_lineage.values[m].sum())}
        for meth, key in (("pLIN v4.1 L5", "pLIN_L5"), ("MOB-suite secondary", "MOB_secondary")):
            pred = same_group(lab[meth], a, b)
            p, r, f = metrics(pred[m], pairs.same_lineage.values[m], w[m])
            row[key] = round(float(f), 3)
        out[name] = row
    json.dump(out, open(os.path.join(C, "taxon_breakdown.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
