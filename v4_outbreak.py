#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Outbreak validation of pLIN v4 (secondary analysis; not pre-registered).

The 74 published outbreak plasmids (27 studies; outbreak_validation/, metadata
in output/outbreak_validation_founder_results.tsv) are typed with
`plin_v4.py query`, exactly as a user would (whole-plasmid match for the 56
already in the database, the protein path for the rest). For every grouping
we report the share of same-study pairs placed together (sensitivity for
outbreak links) and of different-study pairs placed together (false links),
alongside pLIN v3 and MOB-suite.

Usage:
  python v4_outbreak.py
Output: output/backbone_v4/evaluation/outbreak_{codes,metrics}.tsv
"""

import glob
import os
import subprocess
import sys
from itertools import combinations

import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V4 = os.path.join(BASE_DIR, "output", "backbone_v4")
OB_DIR = os.path.join(BASE_DIR, "outbreak_validation")
MOB_DIR = os.path.join(BASE_DIR, "output", "comparator_benchmark", "mob_typer")


def prefix(code, k):
    return ".".join(str(code).split(".")[:k])


def main():
    ob = pd.read_csv(os.path.join(BASE_DIR, "output", "outbreak_validation_founder_results.tsv"), sep="\t")
    files = {os.path.basename(f).rsplit(".", 1)[0]: f
             for f in glob.glob(os.path.join(OB_DIR, "*.fasta")) + glob.glob(os.path.join(OB_DIR, "expanded_sequences", "*.fasta"))}
    work = os.path.join(V4, "outbreak")
    os.makedirs(work, exist_ok=True)
    lst = os.path.join(work, "inputs.txt")
    open(lst, "w").write("\n".join(files[a] for a in ob.accession) + "\n")
    out = os.path.join(work, "v4_query_codes.tsv")
    subprocess.run([sys.executable, os.path.join(BASE_DIR, "plin_v4.py"), "query", lst, "--out", out,
                    "--workdir", os.path.join(work, "query_work"), "--threads", "8"], check=True)
    v4 = pd.read_csv(out, sep="\t").set_index("plasmid_id")
    ob["pLIN_v4"] = ob.accession.map(v4.pLIN_v4)
    ob["v4_matched_database"] = ob.accession.map(v4.matched).notna()
    ob["v4_provisional_from"] = ob.accession.map(v4.get("provisional_from"))

    def mob(acc, col):
        f = os.path.join(MOB_DIR, f"{acc}.txt")
        if not os.path.exists(f):
            return None
        x = pd.read_csv(f, sep="\t", dtype=str).iloc[0].get(col)
        return None if pd.isna(x) or x in ("-", "") else x
    ob["mob_primary"] = [mob(a, "primary_cluster_id") for a in ob.accession]
    ob["mob_secondary"] = [mob(a, "secondary_cluster_id") for a in ob.accession]
    ob.to_csv(os.path.join(V4, "evaluation", "outbreak_codes.tsv"), sep="\t", index=False)

    groupings = {f"pLIN v4 L{k}": ob.pLIN_v4.map(lambda c, k=k: prefix(c, k)) for k in range(1, 7)}
    groupings.update({f"pLIN v3 L{k}": ob.pLIN.map(lambda c, k=k: prefix(c, k)) for k in (3, 5, 6)})
    groupings["MOB-suite primary"] = ob.mob_primary
    groupings["MOB-suite secondary"] = ob.mob_secondary
    groupings["pLIN v4 L6 AND MOB-suite secondary"] = pd.Series(
        [f"{a}|{b}" if b is not None else None for a, b in zip(groupings["pLIN v4 L6"], ob.mob_secondary)])
    rows = []
    for name, g in groupings.items():
        g = g.reset_index(drop=True)
        same, diff = [], []
        for i, j in combinations(range(len(ob)), 2):
            together = g.iloc[i] is not None and pd.notna(g.iloc[i]) and g.iloc[i] == g.iloc[j]
            (same if ob.study.iloc[i] == ob.study.iloc[j] else diff).append(together)
        rows.append({"method": name, "assigned_pct": round(100 * g.notna().mean(), 1),
                     "same_study_pairs": len(same), "same_study_together_pct": round(100 * np.mean(same), 1),
                     "different_study_pairs": len(diff), "different_study_together_pct": round(100 * np.mean(diff), 2),
                     "false_links": int(np.sum(diff))})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(V4, "evaluation", "outbreak_metrics.tsv"), sep="\t", index=False)
    print(f"{ob.v4_matched_database.sum()} of {len(ob)} matched a database plasmid exactly; "
          f"provisional codes: {ob.v4_provisional_from.notna().sum()}")
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
