#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Head-to-head benchmark of pLIN against MOB-suite clusters and pling.

Truth set 1 (alignment): the plasmid pairs aligned by validate_alignment_backbone.py
(single-linkage and founder runs pooled). A pair is
  same lineage      if bidirectional alignment fraction >= 0.8 and identity >= 99%
  related backbone  if >= 50% of the shorter plasmid's backbone aligns.
For each method and grouping level, a pair is "predicted related" when both
plasmids fall in the same group. We report precision, recall and F1, and the
share of plasmids the method could assign to any group. Caveat: pairs were
sampled by pLIN code sharing, so the pair set is enriched for pLIN-related
pairs; all methods are scored on the same pairs.

Truth set 2 (outbreaks): 74 published outbreak plasmids; share of same-study
and different-study pairs placed in the same group.

Methods
  pLIN L6 / L5 / L3            founder codes (output/pLIN_assignments.tsv; outbreak codes
                               from output/outbreak_validation_founder_results.tsv)
  pLIN L6 + protein flag       L6, or backbone-protein containment >= 0.9 with cosine <= 0.005
                               (output/protein_pilot/pair_protein_metrics.tsv)
  MOB-suite primary/secondary  mob_typer cluster IDs (MOB-suite reference database)
  pling subcommunity/community pling 3 cluster align (sourmash prefilter)

Usage:
  python benchmark_vs_comparators.py
Output: output/comparator_benchmark/{benchmark_alignment_truth.tsv, benchmark_outbreak.tsv, summary.json}
"""

import os
import glob
import json
from itertools import combinations

import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CB = os.path.join(BASE_DIR, "output", "comparator_benchmark")
PAIRS = [os.path.join(BASE_DIR, "output", "alignment_validation_single_linkage_v3.2.3", "pair_metrics.tsv"),
         os.path.join(BASE_DIR, "output", "alignment_validation", "pair_metrics.tsv")]


def load_mob():
    rows = []
    for f in glob.glob(os.path.join(CB, "mob_typer", "*.txt")):
        d = pd.read_csv(f, sep="\t", dtype=str)
        if len(d):
            r = d.iloc[0]
            rows.append({"file_id": os.path.basename(f)[:-4], "mob_primary": r.get("primary_cluster_id"),
                         "mob_secondary": r.get("secondary_cluster_id")})
    m = pd.DataFrame(rows).set_index("file_id")
    for c in ("mob_primary", "mob_secondary"):
        m[c] = m[c].where(~m[c].isin(["-", "", None]) & m[c].notna())
    return m


def load_pling(out_dir):
    """plasmid -> community / subcommunity from pling's plasnet outputs.

    Hub plasmids are left out of typing.tsv by pling and so count as unassigned.
    Returns {} if pling has not finished; raises if it finished without typing.
    """
    typing = glob.glob(os.path.join(out_dir, "dcj_thresh_*_graph", "objects", "typing.tsv"))
    if not typing:
        if glob.glob(os.path.join(out_dir, "*_distances.tsv")):
            raise FileNotFoundError(f"pling distances exist but no typing.tsv under {out_dir}")
        return {}
    files = {"community": os.path.join(out_dir, "containment", "containment_communities", "objects", "communities.tsv"),
             "subcommunity": typing[0]}
    res = {}
    for col, path in files.items():
        d = pd.read_csv(path, sep="\t", dtype=str)
        res[col] = dict(zip(d["plasmid"], d[d.columns[1]]))
    return res


def prefix(code, k):
    return ".".join(str(code).split(".")[:k])


def score(same_group, truth):
    tp = int((same_group & truth).sum())
    pred, pos = int(same_group.sum()), int(truth.sum())
    p = tp / pred if pred else np.nan
    r = tp / pos if pos else np.nan
    f1 = 2 * p * r / (p + r) if pred and pos and (p + r) else np.nan
    return {"precision": round(p, 4), "recall": round(r, 4), "F1": round(f1, 4), "TP": tp,
            "predicted_pairs": pred, "true_pairs": pos}


def main():
    pairs = pd.concat([pd.read_csv(p, sep="\t") for p in PAIRS]).drop_duplicates(["plasmid_shorter", "plasmid_longer"])
    pairs = pairs.reset_index(drop=True)
    founder = pd.read_csv(os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv"), sep="\t") \
        .drop_duplicates("plasmid_id").set_index("plasmid_id").pLIN
    prot = pd.read_csv(os.path.join(BASE_DIR, "output", "protein_pilot", "pair_protein_metrics.tsv"), sep="\t") \
        .set_index(["plasmid_shorter", "plasmid_longer"])
    mob = load_mob()
    pling = load_pling(os.path.join(CB, "pling_out"))

    a, b = pairs.plasmid_shorter, pairs.plasmid_longer
    truth = {"same_lineage": (pairs.AF_min >= 0.8) & (pairs.ANI_aln >= 99),
             "related_backbone": pairs.bb_AF_shorter.fillna(pairs.AF_shorter) >= 0.5}

    def same(lookup):
        ga, gb = a.map(lookup), b.map(lookup)
        return (ga.notna() & gb.notna() & (ga == gb)).values, float(pd.concat([ga, gb]).notna().mean())

    methods = {}
    for k, name in ((6, "pLIN L6 (lineage)"), (5, "pLIN L5 (clone group)"), (4, "pLIN L4"), (3, "pLIN L3")):
        methods[name] = same(founder.map(lambda c, k=k: prefix(c, k)).to_dict())
    pc = prot.reindex(list(zip(a, b)))
    flag = (pc.bb_containment.values >= 0.9) & (pairs.cosine_distance.values <= 0.005)
    l6, _ = methods["pLIN L6 (lineage)"]
    methods["pLIN L6 + backbone-protein flag"] = (l6 | flag, 1.0)
    methods["MOB-suite secondary cluster"] = same(mob.mob_secondary.to_dict())
    methods["MOB-suite primary cluster"] = same(mob.mob_primary.to_dict())
    if "subcommunity" in pling:
        methods["pling subcommunity"] = same(pling["subcommunity"])
    if "community" in pling:
        methods["pling community"] = same(pling["community"])

    rows = []
    for name, (sg, cov) in methods.items():
        for tname, t in truth.items():
            rows.append({"method": name, "truth": tname, "plasmids_assigned_pct": round(100 * cov, 1),
                         **score(pd.Series(sg), pd.Series(t.values))})
    al = pd.DataFrame(rows)
    al.to_csv(os.path.join(CB, "benchmark_alignment_truth.tsv"), sep="\t", index=False)
    print(al.to_string(index=False))

    # outbreak set
    ob = pd.read_csv(os.path.join(BASE_DIR, "output", "outbreak_validation_founder_results.tsv"), sep="\t")
    ob["mob_primary"] = ob.accession.map(mob.mob_primary)
    ob["mob_secondary"] = ob.accession.map(mob.mob_secondary)
    groupings = {"pLIN L6": ob.pLIN, "pLIN L5": ob.pLIN.map(lambda c: prefix(c, 5)),
                 "pLIN L3": ob.pLIN.map(lambda c: prefix(c, 3)),
                 "MOB-suite secondary": ob.mob_secondary, "MOB-suite primary": ob.mob_primary}
    orow = []
    for name, g in groupings.items():
        same_s, diff_s = [], []
        for i, j in combinations(range(len(ob)), 2):
            together = pd.notna(g.iloc[i]) and g.iloc[i] == g.iloc[j]
            (same_s if ob.study.iloc[i] == ob.study.iloc[j] else diff_s).append(together)
        orow.append({"method": name, "assigned_pct": round(100 * g.notna().mean(), 1),
                     "same_study_pairs_together_pct": round(100 * np.mean(same_s), 1),
                     "different_study_pairs_together_pct": round(100 * np.mean(diff_s), 2)})
    obdf = pd.DataFrame(orow)
    obdf.to_csv(os.path.join(CB, "benchmark_outbreak.tsv"), sep="\t", index=False)
    print(obdf.to_string(index=False))
    json.dump({"n_pairs": int(len(pairs)), "n_same_lineage": int(truth["same_lineage"].sum()),
               "n_related_backbone": int(truth["related_backbone"].sum()),
               "pling_available": bool(pling)}, open(os.path.join(CB, "summary.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
