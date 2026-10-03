#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Agreement of pLIN v4.1 levels with published plasmid taxonomic units (PTUs;
Redondo-Salvo et al. 2020), alongside MOB-suite, using the scoring of v4_ptu.py
(adjusted Rand index; pairwise precision, recall, F1 for "same PTU").

  confirmatory test set   v4.1 test plasmids that have a PTU (small; descriptive)
  all with a PTU          every release plasmid with a PTU (exploratory)

Usage:
  python v41_ptu.py
Output: output/backbone_v41/ptu_agreement.tsv
"""

import os

import pandas as pd

from v4_ptu import mob_labels, release_ids_with_ptu, score

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V41 = os.path.join(BASE_DIR, "output", "backbone_v41")


def main():
    truth_by_id = release_ids_with_ptu()
    truth = {i.replace("RefSeq_", "", 1): p for i, p in truth_by_id.items()}
    codes = pd.read_csv(os.path.join(V41, "codes_v41.tsv"), sep="\t").set_index("accession").pLIN_v41
    id_of = pd.read_csv(os.path.join(V41, "codes_v41.tsv"), sep="\t").set_index("accession").plasmid_id
    test = set(pd.read_csv(os.path.join(V41, "confirm", "test_plasmids.tsv"), sep="\t").plasmid_id.str.replace("^RefSeq_", "", regex=True))
    sets = {"confirmatory test set": sorted(set(truth) & test & set(codes.index)),
            "all release plasmids with a PTU": sorted(set(truth) & set(codes.index))}
    rows = []
    for name, accs in sets.items():
        lab = {f"pLIN v4.1 L{k}": {a: ".".join(codes[a].split(".")[:k]) for a in accs} for k in range(1, 7)}
        mob = mob_labels([id_of[a] for a in accs])
        mob = {a: mob[id_of[a]] for a in accs if id_of[a] in mob}
        lab["MOB-suite primary"] = {a: m[0] for a, m in mob.items()}
        lab["MOB-suite secondary"] = {a: m[1] for a, m in mob.items()}
        for m, l in lab.items():
            sub = [a for a in accs if not m.startswith("MOB") or a in l]
            rows.append({"set": name, "method": m, **score(l, truth, sub)})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(V41, "ptu_agreement.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
