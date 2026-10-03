#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Side-by-side comparison of v4.1 candidate designs (development data only).

Reads output/backbone_v41/dev/eval_*.json (v41_dev_eval.py) and prints, per
design: the best related-backbone F1 and same-lineage F1 over levels (and at
which level), re-query reproduction, mutant code retention at L3 and L6, the
Swiss VIM-1 checks, outbreak sensitivity / false links at L6, cluster counts.

Usage:
  python v41_dev_compare.py
Output: output/backbone_v41/dev/compare.tsv
"""

import glob
import json
import os

import pandas as pd

DEV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "backbone_v41", "dev")


def main():
    rows = []
    for f in sorted(glob.glob(os.path.join(DEV, "eval_*.json"))):
        r = json.load(open(f))
        name = os.path.basename(f)[5:-5]
        pa = r["pairs"]["all"]
        bb = {k: v[0] for k, v in pa["related_backbone"].items()}
        sl = {k: v[0] for k, v in pa["same_lineage"].items()}
        bbk, slk = max(bb, key=bb.get), max(sl, key=sl.get)
        lv = r["config"]["levels"]
        rows.append({
            "design": name, "engine": r["config"].get("engine", "founder"),
            "levels": " ".join(f"{k}{t}" for k, t in lv),
            "backbone_F1_best": f"{bb[bbk]:.3f}@{bbk}", "lineage_F1_best": f"{sl[slk]:.3f}@{slk}",
            "lineage_F1_L6": sl["L6"], "lineage_P_L6": pa["same_lineage"]["L6"][1], "lineage_R_L6": pa["same_lineage"]["L6"][2],
            "requery": r["requery_identical"],
            "mut0.1%_L3": r["mutants"]["0.001"][2], "mut0.1%_L6": r["mutants"]["0.001"][5],
            "mut0.01%_L6": r["mutants"]["0.0001"][5],
            "vim_samelin_allL": f"{r['vim']['same_lineage_pairs_sharing_all_levels']}/{r['vim']['same_lineage_pairs']}",
            "vim_other_split": f"{r['vim']['other_pairs_split_at_last_level']}/{r['vim']['other_pairs']}",
            "vim_min_shared": r["vim"]["all_vim_pairs_min_shared_level"],
            "outbreak_L6_same%": r["outbreak"]["L6"][0], "outbreak_L6_diff%": r["outbreak"]["L6"][1],
            "clusters_L1": r["clusters"]["L1"], "clusters_L6": r["clusters"]["L6"],
            "build_s": r["build_seconds"]})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DEV, "compare.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 300)
    pd.set_option("display.max_columns", 40)
    print(df.drop(columns=["levels"]).to_string(index=False))
    print("\n" + df[["design", "levels"]].to_string(index=False))


if __name__ == "__main__":
    main()
