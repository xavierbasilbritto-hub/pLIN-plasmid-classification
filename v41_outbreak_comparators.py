#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Comparative outbreak study: pLIN v4.1 vs MOB-suite, pling and mge-cluster on the
74 published outbreak plasmids (27 studies). Development data for pLIN (used during
the v4.1 design); comparators run with default settings.

For each method and level: share of same-study pairs placed together (outbreak
links found) and of different-study pairs placed together (links between outbreaks;
some are real, as several studies describe the same internationally spread
plasmid). Unassigned plasmids are never grouped.

Usage:
  python v41_outbreak_comparators.py
Output: output/backbone_v41/case_studies/outbreak_comparators.tsv
"""

import glob
import os
import subprocess
from itertools import combinations

import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "output", "backbone_v41", "case_studies")
WORK = os.path.join(OUT, "outbreak_comparators")
PLING = os.path.expanduser("~/miniconda3/envs/pling_env/bin")
MGE = os.path.expanduser("~/miniconda3/envs/mge_cluster_env/bin/mge_cluster")
MOB_DIR = os.path.join(BASE_DIR, "output", "comparator_benchmark", "mob_typer")


def run_tools(paths):
    os.makedirs(WORK, exist_ok=True)
    lst = os.path.join(WORK, "inputs.txt")
    open(lst, "w").write("\n".join(paths) + "\n")
    env = dict(os.environ, PATH=PLING + os.pathsep + os.environ["PATH"])
    if not glob.glob(os.path.join(WORK, "pling", "dcj_thresh_*_graph", "objects", "typing.tsv")):
        subprocess.run(["pling", "cluster", "align", "inputs.txt", "pling", "--sourmash", "--cores", "8",
                        "--visualisation", "none"], cwd=WORK, env=env, capture_output=True)
    if not os.path.exists(os.path.join(WORK, "mge", "mge-cluster_results.csv")):
        subprocess.run([MGE, "--create", "--input", lst, "--outdir", os.path.join(WORK, "mge"), "--threads", "8"],
                       capture_output=True)


def main():
    ob = pd.read_csv(os.path.join(BASE_DIR, "output", "outbreak_validation_founder_results.tsv"), sep="\t")
    files = {os.path.basename(f).rsplit(".", 1)[0]: f for f in
             glob.glob(os.path.join(BASE_DIR, "outbreak_validation", "*.fasta")) +
             glob.glob(os.path.join(BASE_DIR, "outbreak_validation", "expanded_sequences", "*.fasta"))}
    run_tools([files[a] for a in ob.accession])
    lab = {}
    v41 = pd.read_csv(os.path.join(OUT, "outbreak_codes.tsv"), sep="\t").set_index("plasmid_id").pLIN_v41
    for k in (1, 3, 5, 6):
        lab[f"pLIN v4.1 L{k}"] = {a: ".".join(v41[a].split(".")[:k]) for a in ob.accession}
    clean = lambda x: None if pd.isna(x) or x in ("-", "") else x
    mob = {a: pd.read_csv(os.path.join(MOB_DIR, f"{a}.txt"), sep="\t", dtype=str).iloc[0] for a in ob.accession}
    lab["MOB-suite primary"] = {a: clean(mob[a].get("primary_cluster_id")) for a in ob.accession}
    lab["MOB-suite secondary"] = {a: clean(mob[a].get("secondary_cluster_id")) for a in ob.accession}
    typing = glob.glob(os.path.join(WORK, "pling", "dcj_thresh_*_graph", "objects", "typing.tsv"))
    if typing:
        t = dict(pd.read_csv(typing[0], sep="\t", dtype=str).iloc[:, :2].values)
        c = dict(pd.read_csv(os.path.join(WORK, "pling", "containment", "containment_communities", "objects",
                                          "communities.tsv"), sep="\t", dtype=str).iloc[:, :2].values)
        lab["pling subcommunity"] = {a: t.get(a) for a in ob.accession}
        lab["pling community"] = {a: c.get(a) for a in ob.accession}
    mres = os.path.join(WORK, "mge", "mge-cluster_results.csv")
    if os.path.exists(mres):
        m = pd.read_csv(mres, dtype=str)
        mg = {s: (None if c.strip() == "-1" else c.strip()) for s, c in zip(m.Sample_Name, m.Standard_Cluster)}
        lab["mge-cluster"] = {a: mg.get(a) for a in ob.accession}
    st = dict(zip(ob.accession, ob.study))
    rows = []
    for name, l in lab.items():
        same, diff = [], []
        for a, b in combinations(ob.accession, 2):
            tog = l[a] is not None and l[a] == l[b]
            (same if st[a] == st[b] else diff).append(tog)
        rows.append({"method": name, "assigned_pct": round(100 * np.mean([l[a] is not None for a in ob.accession]), 1),
                     "same_study_pairs": len(same), "same_study_together_pct": round(100 * np.mean(same), 1),
                     "different_study_pairs": len(diff), "different_study_together_pct": round(100 * np.mean(diff), 2),
                     "links_between_outbreaks": int(np.sum(diff))})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "outbreak_comparators.tsv"), sep="\t", index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
