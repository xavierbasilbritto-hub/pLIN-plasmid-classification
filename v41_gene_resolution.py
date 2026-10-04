#!/usr/bin/env python3
# Copyright (C) 2025-2026 Basil Britto Xavier. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Do pLIN codes predict the resistance gene better than replicon (Inc/Rep) types?

Replicon typing is the standard way to name a resistance plasmid, but a replicon type is a
property of the replication machinery, not of the cargo: one Inc type carries many different
resistance genes. If a pLIN lineage is a far better predictor of which gene a plasmid carries,
that is a practical argument for the nomenclature that replicon typing cannot answer.

Two measures on release plasmids matched to PLSDB AMRFinderPlus calls:
  purity  for every group of at least 10 carriers, the share carrying the single commonest key
          gene (carbapenemase, mcr, CTX-M or tet(X), as in v4_onehealth.py)
  AMI     adjusted mutual information between the grouping and the resistance gene, over
          plasmids carrying exactly one key gene, so the label is unambiguous

Usage:
  python v41_gene_resolution.py
Output: output/backbone_v41/onehealth/gene_resolution.json and gene_resolution_examples.tsv
"""

import json
import os

import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score

from v4_onehealth import AMR_GROUPS, META

BASE = os.path.dirname(os.path.abspath(__file__))
V41 = os.path.join(BASE, "output", "backbone_v41")
OH = os.path.join(V41, "onehealth")
MIN_CARRIERS = 10


def load():
    sec = pd.read_csv(os.path.join(OH, "plasmid_sectors.tsv"), sep="\t", dtype=str)
    rel = pd.read_csv(os.path.join(V41, "release", "plin_v41_codes.tsv.gz"), sep="\t",
                      usecols=["accession", "pLIN_v41", "inc_type"])
    amr = pd.read_csv(os.path.join(META, "amr.tsv"), sep="\t",
                      usecols=["NUCCORE_ACC", "gene_symbol"], dtype=str)
    key = pd.concat([amr[amr.gene_symbol.str.contains(rx, regex=True, na=False)].assign(group=g)
                     for g, rx in AMR_GROUPS.items()]).drop_duplicates(["NUCCORE_ACC", "gene_symbol"])
    m = sec[["NUCCORE_ACC", "accession"]].merge(rel, on="accession")
    carriers = key.merge(m, on="NUCCORE_ACC")
    for df in (m, carriers):
        df["L1"] = df.pLIN_v41.str.split(".").str[0]
        df["L5"] = df.pLIN_v41.map(lambda c: ".".join(str(c).split(".")[:5]))
    return m, key, carriers


def purity(carriers, col):
    """Share of each group's carriers that have the group's commonest key gene."""
    per = carriers.groupby([col, "gene_symbol"]).NUCCORE_ACC.nunique().reset_index()
    total = carriers.groupby(col).NUCCORE_ACC.nunique()
    big = total[total >= MIN_CARRIERS].index
    per = per[per[col].isin(big)]
    top = per.sort_values("NUCCORE_ACC").groupby(col).tail(1).set_index(col)
    share = top.NUCCORE_ACC / total[top.index]
    return {"groups": int(len(share)), "median_dominant_pct": round(100 * float(share.median()), 1),
            "pct_groups_at_least_75": round(100 * float((share >= 0.75).mean()), 1)}


def main():
    m, key, carriers = load()
    # unambiguous subset: plasmids carrying exactly one key gene
    n_genes = key.groupby("NUCCORE_ACC").gene_symbol.nunique()
    single = key[key.NUCCORE_ACC.isin(n_genes[n_genes == 1].index)].merge(m, on="NUCCORE_ACC")

    levels = [("inc_type", "Inc/Rep type"), ("L1", "pLIN L1 family"),
              ("L5", "pLIN L5 lineage"), ("pLIN_v41", "pLIN L6 clone")]
    out = {"plasmids_with_metadata": int(len(m)), "carriers": int(carriers.NUCCORE_ACC.nunique()),
           "single_gene_plasmids": int(len(single)), "distinct_genes": int(single.gene_symbol.nunique()),
           "min_carriers_per_group": MIN_CARRIERS, "levels": {}}
    for col, lab in levels:
        s = single.dropna(subset=[col])
        out["levels"][lab] = {**purity(carriers.dropna(subset=[col]), col),
                              "AMI_with_resistance_gene": round(
                                  float(adjusted_mutual_info_score(s[col], s.gene_symbol)), 3)}

    # worked example: one replicon type, many genes; its pLIN lineages, one gene each
    inc = carriers.groupby("inc_type").NUCCORE_ACC.nunique().idxmax()
    sub = carriers[carriers.inc_type == inc]
    rows = []
    total = sub.groupby("L5").NUCCORE_ACC.nunique()
    for code, n in total[total >= 20].sort_values(ascending=False).head(8).items():
        g = sub[sub.L5 == code].gene_symbol.value_counts()
        rows.append({"inc_type": inc, "pLIN_L5": code, "plasmids": int(n),
                     "dominant_gene": g.index[0], "dominant_pct": round(100 * g.iloc[0] / n, 1)})
    ex = pd.DataFrame(rows)
    ex.to_csv(os.path.join(OH, "gene_resolution_examples.tsv"), sep="\t", index=False)
    out["example"] = {"inc_type": inc, "carriers": int(sub.NUCCORE_ACC.nunique()),
                      "distinct_key_genes": int(sub.gene_symbol.nunique()),
                      "lineages_shown": len(ex)}
    json.dump(out, open(os.path.join(OH, "gene_resolution.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))
    print(ex.to_string(index=False))


if __name__ == "__main__":
    main()
