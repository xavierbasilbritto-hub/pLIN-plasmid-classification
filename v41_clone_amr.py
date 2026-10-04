#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Key resistance genes in the largest pLIN v4.1 outbreak clones (L6), descriptive.

Release plasmids matched to PLSDB 2024_05_31_v2 metadata carry PLSDB's AMRFinderPlus calls. For the
20 largest L6 clones that carry at least one key resistance gene (carbapenemase, mcr, CTX-M or
tet(X); same groups as v4_onehealth.py), the share of each clone's plasmids that carry each key gene.

Usage:
  python v41_clone_amr.py
Output: output/backbone_v41/onehealth/clone_amr_matrix.tsv (rows = clones, columns = genes, % carrying),
        output/backbone_v41/onehealth/clone_amr_summary.json
"""

import json
import os

import pandas as pd

from v4_onehealth import AMR_GROUPS, META

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OH = os.path.join(BASE_DIR, "output", "backbone_v41", "onehealth")
N_CLONES, MIN_GENE_PCT = 20, 20


def main():
    sec = pd.read_csv(os.path.join(OH, "plasmid_sectors.tsv"), sep="\t", dtype=str)
    sec = sec.rename(columns={"pLIN_v4": "pLIN_v41"})             # column holds the v4.1 codes (see build_facts.py)
    sec["L6"] = sec.pLIN_v41
    a = pd.read_csv(os.path.join(META, "amr.tsv"), sep="\t", usecols=["NUCCORE_ACC", "gene_symbol"], dtype=str)
    key = pd.concat([a[a.gene_symbol.str.contains(rx, regex=True, na=False)].assign(group=g)
                     for g, rx in AMR_GROUPS.items()]).drop_duplicates(["NUCCORE_ACC", "gene_symbol"])
    size = sec.groupby("L6").size()
    carriers = key.merge(sec[["NUCCORE_ACC", "L6"]], on="NUCCORE_ACC")
    with_key = carriers.L6.unique()
    top = size[size.index.isin(with_key)].sort_values(ascending=False).head(N_CLONES)
    sub = carriers[carriers.L6.isin(top.index)]
    pct = (sub.groupby(["L6", "gene_symbol"]).NUCCORE_ACC.nunique().unstack(fill_value=0)
           .div(top, axis=0) * 100).loc[top.index]
    keep = pct.columns[(pct >= MIN_GENE_PCT).any()]
    pct = pct[keep].round(1)
    pct.insert(0, "plasmids", top)
    pct.index.name = "pLIN_L6_clone"
    pct.to_csv(os.path.join(OH, "clone_amr_matrix.tsv"), sep="\t")
    summary = {"plasmids_with_metadata": int(len(sec)), "L6_clones_with_key_gene": int(len(with_key)),
               "clones_shown": int(len(top)), "genes_shown": int(len(keep)), "min_gene_pct": MIN_GENE_PCT,
               "largest_clone_plasmids": int(top.iloc[0])}
    json.dump(summary, open(os.path.join(OH, "clone_amr_summary.json"), "w"), indent=1)
    print(pct.to_string())
    print(summary)


if __name__ == "__main__":
    main()
