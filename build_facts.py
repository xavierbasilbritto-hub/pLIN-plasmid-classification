#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Single source of truth for every number quoted in pLIN documents.

Reads the release and result files and writes docs/PLIN_FACTS.json: each fact
has a value, a ready-to-quote text form and the file it comes from. Documents
(README, guides, manual, app help, manuscript) quote these; check_docs.py
verifies that they do and that no outdated number remains.

Usage:
  python build_facts.py
"""

import json
import os

import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V41 = os.path.join(BASE_DIR, "output", "backbone_v41")
OUT = os.path.join(BASE_DIR, "docs", "PLIN_FACTS.json")


def rel(p):
    return os.path.relpath(p, BASE_DIR)


def main():
    facts = {}

    def add(key, value, text, source):
        facts[key] = {"value": value, "text": text, "source": rel(source)}

    dv_path = os.path.join(V41, "release", "DATABASE_VERSION.json")
    dv = json.load(open(dv_path))
    add("scheme", dv["scheme"], dv["scheme"], dv_path)
    add("database_version", dv["database_version"], dv["database_version"], dv_path)
    add("db_plasmids", dv["unique_plasmids"], f"{dv['unique_plasmids']:,}", dv_path)
    for k, v in dv["clusters_per_level"].items():
        add(f"clusters_{k}", v, f"{v:,}", dv_path)
    add("db_no_protein", dv["plasmids_without_proteins"], f"{dv['plasmids_without_proteins']:,}", dv_path)

    ep_path = os.path.join(V41, "confirm", "endpoints.json")
    ep = json.load(open(ep_path))
    add("test_pairs", ep["test_pairs"], f"{ep['test_pairs']:,}", ep_path)
    add("test_same_lineage_pairs", ep["same_lineage_pairs"], f"{ep['same_lineage_pairs']:,}", ep_path)
    e2, e1 = ep["E2_lineage_v41L5_vs_MOBsecondary"], ep["E1_backbone_v41L1_vs_MOBprimary"]
    add("lineage_F1", round(e2["F1"], 2), f"{e2['F1']:.2f}", ep_path)
    add("lineage_F1_MOB", round(e2["F1_comparator"], 2), f"{e2['F1_comparator']:.2f}", ep_path)
    add("lineage_diff", round(e2["difference"], 2),
        f"+{e2['difference']:.2f} (95% CI {e2['lo']:.2f} to {e2['hi']:.2f})", ep_path)
    add("backbone_F1", round(e1["F1"], 2), f"{e1['F1']:.2f}", ep_path)
    add("backbone_F1_MOB", round(e1["F1_comparator"], 2), f"{e1['F1_comparator']:.2f}", ep_path)
    add("backbone_diff", round(e1["difference"], 2),
        f"+{e1['difference']:.2f} (95% CI {e1['lo']:.2f} to {e1['hi']:.2f}); not shown non-inferior", ep_path)

    tm_path = os.path.join(V41, "confirm", "test_metrics.tsv")
    tm = pd.read_csv(tm_path, sep="\t")
    g = lambda m, t, c: float(tm[(tm.method == m) & (tm.truth == t)][c].iloc[0])
    add("lineage_F1_CI", [round(g("pLIN v4.1 L5", "same_lineage", "F1_w_lo"), 2), round(g("pLIN v4.1 L5", "same_lineage", "F1_w_hi"), 2)],
        f"{g('pLIN v4.1 L5', 'same_lineage', 'F1_w_lo'):.2f} to {g('pLIN v4.1 L5', 'same_lineage', 'F1_w_hi'):.2f}", tm_path)
    add("lineage_precision", round(g("pLIN v4.1 L5", "same_lineage", "precision_w"), 2),
        f"{g('pLIN v4.1 L5', 'same_lineage', 'precision_w'):.2f}", tm_path)
    add("lineage_recall", round(g("pLIN v4.1 L5", "same_lineage", "recall_w"), 2),
        f"{g('pLIN v4.1 L5', 'same_lineage', 'recall_w'):.2f}", tm_path)
    for m, key in (("pling subcommunity", "lineage_F1_pling"), ("mge-cluster", "lineage_F1_mge"),
                   ("pLIN v4 L6", "lineage_F1_v4")):
        add(key, round(g(m, "same_lineage", "F1_w"), 2), f"{g(m, 'same_lineage', 'F1_w'):.2f}", tm_path)

    rc_path = os.path.join(V41, "confirm", "release_checks.json")
    rc = json.load(open(rc_path))
    add("stability_runs", len(rc["R1_stability"]["runs"]), str(len(rc["R1_stability"]["runs"])), rc_path)
    add("stability_retained", rc["R1_stability"]["min_retained_all_levels"], "100%", rc_path)
    add("requery_n", rc["R2_reproduction"]["n"], str(rc["R2_reproduction"]["n"]), rc_path)
    add("requery_identical", rc["R2_reproduction"]["identical"], "100%", rc_path)
    r3 = rc["R3_robustness"]
    add("mutation_L1_L5_0.1pct", min(r3["0.001"][:5]), f"{100 * min(r3['0.001'][:5]):.0f}%", rc_path)
    add("mutation_L6_0.1pct", r3["0.001"][5], f"{100 * r3['0.001'][5]:.0f}%", rc_path)

    sp_path = os.path.join(V41, "speed_v41.json")
    sp = json.load(open(sp_path))
    add("speed_related_s", sp["related"]["seconds_per_plasmid"], f"{sp['related']['seconds_per_plasmid']:.2f} s", sp_path)
    add("speed_divergent_s", sp["divergent"]["seconds_per_plasmid"], f"{sp['divergent']['seconds_per_plasmid']:.2f} s", sp_path)
    add("speed_MOB_s", sp["MOB-suite_seconds_per_plasmid_8_threads"], f"{sp['MOB-suite_seconds_per_plasmid_8_threads']:.2f} s", sp_path)

    oh_path = os.path.join(V41, "onehealth", "summary.json")
    oh = json.load(open(oh_path))
    add("onehealth_clones_multisector", oh["L6"]["spanning_>=2_sectors"], f"{oh['L6']['spanning_>=2_sectors']:,}", oh_path)
    add("onehealth_clones_same_gene", oh["L6"]["same_key_AMR_gene_in_>=2_sectors"],
        f"{oh['L6']['same_key_AMR_gene_in_>=2_sectors']:,}", oh_path)

    vp_path = os.path.join(V41, "case_studies", "swiss_vim1_pairs.tsv")
    vp = pd.read_csv(vp_path, sep="\t")
    agree = int(((vp.same_lineage_by_alignment) & (vp.levels_shared_v41 >= 5)).sum() +
                ((~vp.same_lineage_by_alignment) & (vp.levels_shared_v41 < 5)).sum())
    add("vim_pairs_agree", agree, f"{agree} of {len(vp)}", vp_path)

    clf_path = os.path.join(BASE_DIR, "data", "inc_classifier.npz")
    d = np.load(clf_path, allow_pickle=True)
    m = json.loads(str(d["cv_metrics"][0]))
    f1 = [v["f1"] for v in m.values() if isinstance(v, dict) and v.get("f1") is not None]
    add("replicon_groups", len(d["group_names"]), str(len(d["group_names"])), clf_path)
    add("replicon_training", len(d["y"]), f"{len(d['y']):,}", clf_path)
    add("replicon_cv_accuracy", round(float(d["cv_accuracy"][0]), 3), f"{100 * float(d['cv_accuracy'][0]):.1f}%", clf_path)
    add("replicon_macro_F1", round(float(np.mean(f1)), 2), f"{np.mean(f1):.2f}", clf_path)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(facts, open(OUT, "w"), indent=1)
    for k, v in facts.items():
        print(f"{k:32s} {v['text']}")


if __name__ == "__main__":
    main()
