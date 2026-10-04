#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
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
    codes_path = os.path.join(V41, "release", "plin_v41_codes.tsv.gz")
    lens = pd.read_csv(codes_path, sep="\t", usecols=["length_bp"]).length_bp
    add("db_len_min", int(lens.min()), f"{lens.min():,} bp", codes_path)
    add("db_len_max", int(lens.max()), f"{lens.max() / 1e6:.1f} Mb", codes_path)
    from plin_kmers import adaptive_scaled
    add("sketch_small_bp", 5000, "5 kb", "plin_kmers.py")
    add("sketch_small_n", 5000 // adaptive_scaled(5000), str(5000 // adaptive_scaled(5000)), "plin_kmers.py")
    se_path = os.path.join(V41, "size_example.json")
    if os.path.exists(se_path):
        se = json.load(open(se_path))
        add("size_ex_whole", se["source_bp"], f"{se['source_bp'] / 1000:.0f} kb", se_path)
        add("size_ex_segment", se["segment_bp"], f"{se['segment_bp'] // 1000} kb", se_path)
        add("size_ex_cont", se["containment_segment_in_whole"], f"{se['containment_segment_in_whole']:.2f}", se_path)
        add("size_ex_sym", se["symmetric_similarity"], f"{se['symmetric_similarity']:.2f}", se_path)

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
    add("backbone_diff_ci", round(e1["difference"], 2), f"+{e1['difference']:.2f} (95% CI {e1['lo']:.2f} to {e1['hi']:.2f})", ep_path)
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
    add("replicon_IncLM_F1", round(m["IncLM"]["f1"], 2), f"{m['IncLM']['f1']:.2f}", clf_path)

    # does a pLIN lineage predict the resistance gene better than a replicon type? (v41_gene_resolution.py)
    gr_path = os.path.join(V41, "onehealth", "gene_resolution.json")
    if os.path.exists(gr_path):
        gr = json.load(open(gr_path))
        lv = gr["levels"]
        add("gene_res_carriers", gr["carriers"], f"{gr['carriers']:,}", gr_path)
        add("gene_res_single", gr["single_gene_plasmids"], f"{gr['single_gene_plasmids']:,}", gr_path)
        add("gene_res_genes", gr["distinct_genes"], str(gr["distinct_genes"]), gr_path)
        for key, lab in (("inc", "Inc/Rep type"), ("L5", "pLIN L5 lineage"), ("L6", "pLIN L6 clone")):
            add(f"gene_res_{key}_ami", lv[lab]["AMI_with_resistance_gene"],
                f"{lv[lab]['AMI_with_resistance_gene']:.2f}", gr_path)
            add(f"gene_res_{key}_pure", lv[lab]["median_dominant_pct"],
                f"{lv[lab]['median_dominant_pct']:.0f}%", gr_path)
            add(f"gene_res_{key}_groups75", lv[lab]["pct_groups_at_least_75"],
                f"{lv[lab]['pct_groups_at_least_75']:.0f}%", gr_path)
        add("gene_res_example_inc", gr["example"]["inc_type"], gr["example"]["inc_type"], gr_path)
        add("gene_res_example_genes", gr["example"]["distinct_key_genes"],
            str(gr["example"]["distinct_key_genes"]), gr_path)
        add("gene_res_example_n", gr["example"]["carriers"], f"{gr['example']['carriers']:,}", gr_path)

    cs = os.path.join(V41, "case_studies")
    loo_path = os.path.join(cs, "loo_metrics.tsv")
    loo = pd.read_csv(loo_path, sep="\t").set_index("level")
    add("loo_n", int(loo.at["L5", "plasmids_in_database"]), str(int(loo.at["L5", "plasmids_in_database"])), loo_path)
    for k in ("L5", "L6"):
        add(f"loo_recovery_{k}", loo.at[k, "recovery_pct"], f"{loo.at[k, 'recovery_pct']:.1f}%", loo_path)
        add(f"loo_same_pairs_{k}", loo.at[k, "same_study_pairs_pct"], f"{loo.at[k, 'same_study_pairs_pct']:.1f}%", loo_path)
        add(f"loo_diff_pairs_{k}", loo.at[k, "diff_study_pairs_pct"], f"{loo.at[k, 'diff_study_pairs_pct']:.2f}%", loo_path)

    oc_path = os.path.join(cs, "outbreak_comparators.tsv")
    oc = pd.read_csv(oc_path, sep="\t").set_index("method")
    for m, key in (("pLIN v4.1 L5", "ob_L5"), ("pLIN v4.1 L6", "ob_L6"), ("MOB-suite secondary", "ob_MOBsec"),
                   ("MOB-suite primary", "ob_MOBprim"), ("pling subcommunity", "ob_pling_sub"), ("mge-cluster", "ob_mge")):
        add(f"{key}_same", oc.at[m, "same_study_together_pct"], f"{oc.at[m, 'same_study_together_pct']:.1f}%", oc_path)
        add(f"{key}_diff", oc.at[m, "different_study_together_pct"], f"{oc.at[m, 'different_study_together_pct']:.2f}%", oc_path)

    pr_path = os.path.join(cs, "prospective_metrics.tsv")
    pr = pd.read_csv(pr_path, sep="\t").set_index("method")
    add("prosp_n", int(pr.at["pLIN v4.1 L5", "isolates_with_earlier_outbreak_isolate"]),
        str(int(pr.at["pLIN v4.1 L5", "isolates_with_earlier_outbreak_isolate"])), pr_path)
    for m, key in (("pLIN v4.1 L5", "prosp_L5"), ("MOB-suite secondary", "prosp_MOBsec"), ("MOB-suite primary", "prosp_MOBprim")):
        add(f"{key}_warning", pr.at[m, "early_warning_pct"], f"{pr.at[m, 'early_warning_pct']:.1f}%", pr_path)
        add(f"{key}_false_alarm", pr.at[m, "false_alarm_pct"], f"{pr.at[m, 'false_alarm_pct']:.1f}%", pr_path)

    sr_path = os.path.join(V41, "shortread", "fragmentation_metrics.tsv")
    sr = pd.read_csv(sr_path, sep="\t")
    assembled = sr[~sr.group.isin(["all", "not assembled"])]
    add("sr_assembled_n", int(assembled.plasmids.sum()), str(int(assembled.plasmids.sum())), sr_path)
    add("sr_L1_L5_kept", float(assembled.L5_kept_pct.min()), f"{assembled.L5_kept_pct.min():.0f}%", sr_path)
    vb_path = os.path.join(V41, "shortread", "vim_metrics.tsv")
    vb = pd.read_csv(vb_path, sep="\t").set_index("set")
    good = vb.loc[">= 90% of k-mers in one MOB-recon bin"]
    add("sr_vim_reconstructed_n", int(good.plasmids), str(int(good.plasmids)), vb_path)
    add("sr_vim_reconstructed_L5", float(good.L5_kept_pct), f"{good.L5_kept_pct:.0f}%", vb_path)

    dvs_path = os.path.join(V41, "diversity.json")
    dvs = json.load(open(dvs_path))
    add("simpson_replicon", dvs["replicon_typing"], f"{dvs['replicon_typing']:.3f}", dvs_path)
    add("simpson_L5", dvs["pLIN_L5"], f"{dvs['pLIN_L5']:.4f}", dvs_path)
    add("simpson_L6", dvs["pLIN_L6"], f"{dvs['pLIN_L6']:.4f}", dvs_path)

    bs_path = os.path.join(V41, "build_stats.json")
    bs = json.load(open(bs_path))
    for k in ("source_PLSDB_2024_05_31", "source_training_not_in_PLSDB", "source_other_NCBI", "accessions_listed_twice",
              "unique_protein_sequences", "protein_families_in_use", "kmer_hashes_total"):
        add(f"build_{k}", bs[k], f"{bs[k]:,}", bs_path)

    ptu_path = os.path.join(V41, "ptu_agreement.tsv")
    ptu = pd.read_csv(ptu_path, sep="\t")
    allp = ptu[ptu.set == "all release plasmids with a PTU"].set_index("method")
    add("ptu_n", int(allp.at["pLIN v4.1 L1", "n"]), f"{int(allp.at['pLIN v4.1 L1', 'n']):,}", ptu_path)
    add("ptu_ARI_L1", round(allp.at["pLIN v4.1 L1", "ARI"], 2), f"{allp.at['pLIN v4.1 L1', 'ARI']:.2f}", ptu_path)
    add("ptu_ARI_MOBprim", round(allp.at["MOB-suite primary", "ARI"], 2), f"{allp.at['MOB-suite primary', 'ARI']:.2f}", ptu_path)

    # ---- further facts quoted in the manuscript --------------------------------------------
    pct = lambda v, d=0: f"{v:.{d}f}%"
    for m, t, key in (("MOB-suite secondary", "same_lineage", "lineage_MOBsec"), ("pling subcommunity", "related_backbone", "backbone_pling"),
                      ("mge-cluster", "related_backbone", "backbone_mge"), ("pLIN v4 L1", "related_backbone", "backbone_v4"),
                      ("pLIN v4.1 L1", "related_backbone", "backbone_v41L1"), ("MOB-suite primary", "related_backbone", "backbone_MOBprim")):
        for c, s in (("F1_w", "F1"), ("precision_w", "precision"), ("recall_w", "recall"), ("F1_w_lo", "lo"), ("F1_w_hi", "hi")):
            add(f"{key}_{s}", round(g(m, t, c), 2), f"{g(m, t, c):.2f}", tm_path)
    add("test_related_backbone_pairs", ep["related_backbone_pairs"], f"{ep['related_backbone_pairs']:,}", ep_path)
    tp_path = os.path.join(V41, "confirm", "test_plasmids.tsv")
    n_tp = len(pd.read_csv(tp_path, sep="\t"))
    add("test_plasmids", n_tp, f"{n_tp:,}", tp_path)
    for rate, tag in (("0.0001", "0.01pct"), ("0.001", "0.1pct"), ("0.01", "1pct")):
        for k in range(6):
            add(f"mutation_L{k + 1}_{tag}", r3[rate][k], pct(100 * r3[rate][k]), rc_path)
    mp_path = os.path.join(V41, "confirm", "mutation_plasmids.tsv")
    n_mp = len(pd.read_csv(mp_path, sep="\t"))
    add("mutation_n", n_mp, str(n_mp), mp_path)
    add("speed_known_s", sp["known"]["seconds_per_plasmid"], "< 0.01 s", sp_path)
    add("speed_threads", sp["threads"], str(sp["threads"]), sp_path)
    add("speed_n", sp["related"]["plasmids"], str(sp["related"]["plasmids"]), sp_path)

    ps_path = os.path.join(BASE_DIR, "output", "backbone_v4", "evaluation", "stability_pling.tsv")
    ps = pd.read_csv(ps_path, sep="\t").set_index(["mode", "level"])
    add("pling_add_sub_split", ps.at[("add", "subcommunity"), "pairs_split_pct"], pct(ps.at[("add", "subcommunity"), "pairs_split_pct"]), ps_path)
    add("pling_add_n", int(ps.at[("add", "subcommunity"), "n_plasmids"]), f"{int(ps.at[('add', 'subcommunity'), 'n_plasmids']):,}", ps_path)
    ms_path = os.path.join(BASE_DIR, "output", "backbone_v4", "evaluation", "stability_mge_cluster.tsv")
    ms = pd.read_csv(ms_path, sep="\t").set_index("mode")
    add("mge_existing_split", ms.at["existing", "pairs_split_pct"], pct(ms.at["existing", "pairs_split_pct"], 2), ms_path)
    add("mge_rebuild_split", ms.at["rebuild", "pairs_split_pct"], pct(ms.at["rebuild", "pairs_split_pct"]), ms_path)

    ob_path = os.path.join(BASE_DIR, "output", "outbreak_validation_founder_results.tsv")
    ob = pd.read_csv(ob_path, sep="\t")
    add("ob_n", len(ob), str(len(ob)), ob_path)
    add("ob_studies", ob.study.nunique(), str(ob.study.nunique()), ob_path)
    add("ob_regions", ob.country.nunique(), str(ob.country.nunique()), ob_path)
    add("ob_same_pairs", int(oc.at["pLIN v4.1 L5", "same_study_pairs"]), str(int(oc.at["pLIN v4.1 L5", "same_study_pairs"])), oc_path)
    add("ob_diff_pairs", int(oc.at["pLIN v4.1 L5", "different_study_pairs"]), f"{int(oc.at['pLIN v4.1 L5', 'different_study_pairs']):,}", oc_path)
    for m, key in (("pLIN v4.1 L1", "ob_L1"), ("pLIN v4.1 L3", "ob_L3"), ("pling community", "ob_pling_comm")):
        add(f"{key}_same", oc.at[m, "same_study_together_pct"], f"{oc.at[m, 'same_study_together_pct']:.1f}%", oc_path)
        add(f"{key}_diff", oc.at[m, "different_study_together_pct"], f"{oc.at[m, 'different_study_together_pct']:.2f}%", oc_path)
    add("loo_partner_n", int(loo.at["L5", "plasmids_with_outbreak_partner"]), str(int(loo.at["L5", "plasmids_with_outbreak_partner"])), loo_path)
    for k in ("L1", "L3", "L4"):
        add(f"loo_recovery_{k}", loo.at[k, "recovery_pct"], f"{loo.at[k, 'recovery_pct']:.1f}%", loo_path)
    for m, key in (("pLIN v4.1 L1", "prosp_L1"), ("pLIN v4.1 L3", "prosp_L3"), ("pLIN v4.1 L6", "prosp_L6")):
        add(f"{key}_warning", pr.at[m, "early_warning_pct"], f"{pr.at[m, 'early_warning_pct']:.1f}%", pr_path)
        add(f"{key}_false_alarm", pr.at[m, "false_alarm_pct"], f"{pr.at[m, 'false_alarm_pct']:.1f}%", pr_path)
    pa_path = os.path.join(cs, "prospective_arrivals.tsv")
    pa = pd.read_csv(pa_path, sep="\t")
    add("prosp_first_year", int(str(pa.date.min())[:4]), str(pa.date.min())[:4], pa_path)
    add("prosp_last_year", int(str(pa.date.max())[:4]), str(pa.date.max())[:4], pa_path)
    add("prosp_db_first", int(pa.database_plasmids_available.iloc[0]), f"{int(pa.database_plasmids_available.iloc[0]):,}", pa_path)
    add("prosp_db_last", int(pa.database_plasmids_available.iloc[-1]), f"{int(pa.database_plasmids_available.iloc[-1]):,}", pa_path)

    fc_path = os.path.join(V41, "shortread", "fragmentation_codes.tsv")
    fc = pd.read_csv(fc_path, sep="\t")
    un = fc[fc.contigs == 0]
    add("sr_n", len(fc), str(len(fc)), fc_path)
    add("sr_unassembled_n", len(un), str(len(un)), fc_path)
    add("sr_unassembled_range", [int(un.complete_length.min()), int(un.complete_length.max())],
        f"{int(un.complete_length.min())} to {int(un.complete_length.max())} bp", fc_path)
    srg = sr.set_index("group")
    add("sr_L6_1contig", srg.at["1 contig", "L6_kept_pct"], f"{srg.at['1 contig', 'L6_kept_pct']:.1f}%", sr_path)
    add("sr_L6_gt20", srg.at[">20", "L6_kept_pct"], f"{srg.at['>20', 'L6_kept_pct']:.1f}%", sr_path)
    l6a = float((assembled.L6_kept_pct * assembled.plasmids).sum() / assembled.plasmids.sum())
    add("sr_L6_assembled", round(l6a, 1), f"{l6a:.1f}%", sr_path)
    allv = vb.loc["all long-read plasmids"]
    bad = vb.loc["< 90% (binning split or merged the plasmid)"]
    add("sr_vim_n", int(allv.plasmids), str(int(allv.plasmids)), vb_path)
    add("sr_vim_split_n", int(bad.plasmids), str(int(bad.plasmids)), vb_path)
    add("sr_vim_split_L1", float(bad.L1_kept_pct), f"{bad.L1_kept_pct:.1f}%", vb_path)
    add("sr_vim_split_L4", float(bad.L4_kept_pct), f"{bad.L4_kept_pct:.0f}%", vb_path)
    add("sr_vim_reconstructed_L6", float(good.L6_kept_pct), f"{good.L6_kept_pct:.1f}%", vb_path)

    vr_path = os.path.join(V41, "shortread", "vim_real_metrics.tsv")
    vr = pd.read_csv(vr_path, sep="\t").set_index("set")
    rg = vr.loc[">= 90% of k-mers in one MOB-recon bin"]
    rv = vr.loc["blaVIM-1 plasmids"]
    add("sr_real_reconstructed_n", int(rg.plasmids), str(int(rg.plasmids)), vr_path)
    add("sr_real_L1", float(rg.L1_kept_pct), f"{rg.L1_kept_pct:.0f}%", vr_path)
    add("sr_real_L5", float(rg.L5_kept_pct), f"{round(rg.L5_kept_pct / 100 * rg.plasmids):.0f} of {int(rg.plasmids)}", vr_path)
    add("sr_real_L6", float(rg.L6_kept_pct), f"{rg.L6_kept_pct:.1f}%", vr_path)
    add("sr_real_vim_L5", float(rv.L5_kept_pct), f"{round(rv.L5_kept_pct / 100 * rv.plasmids):.0f} of {int(rv.plasmids)}", vr_path)
    rs = vr[vr.index.str.startswith("< 90%")].iloc[0]
    add("sr_real_split_n", int(rs.plasmids), str(int(rs.plasmids)), vr_path)
    add("sr_real_split_L1", float(rs.L1_kept_pct), f"{rs.L1_kept_pct:.0f}%", vr_path)
    add("sr_real_split_L4", float(rs.L4_kept_pct), f"{rs.L4_kept_pct:.0f}%", vr_path)
    vb_bins = pd.read_csv(os.path.join(V41, "shortread", "vim_real_bins.tsv"), sep="\t")
    add("sr_real_vim_isolates", vb_bins.contig_id.str.split("_").str[0].nunique(), str(vb_bins.contig_id.str.split("_").str[0].nunique()), vr_path)

    add("onehealth_annotated", oh["matched_to_PLSDB_2024_05_31"], f"{oh['matched_to_PLSDB_2024_05_31']:,}", oh_path)
    add("onehealth_L6_lineages", oh["L6"]["lineages_with_>=2_annotated_plasmids"], f"{oh['L6']['lineages_with_>=2_annotated_plasmids']:,}", oh_path)
    ex_path = os.path.join(V41, "onehealth", "cross_sector_examples.tsv")
    ex = pd.read_csv(ex_path, sep="\t")
    ex = ex[ex.level == "L6"].sort_values("plasmids_with_gene", ascending=False)
    for gene, key in (("blaOXA-48", "oh_oxa48"), ("blaNDM-5", "oh_ndm5"), ("mcr-1.1", "oh_mcr1")):
        r = ex[ex.gene == gene].iloc[0]
        add(f"{key}_n", int(r.plasmids_with_gene), str(int(r.plasmids_with_gene)), ex_path)
        add(f"{key}_size", int(r.lineage_size), str(int(r.lineage_size)), ex_path)
        add(f"{key}_countries", int(r.countries), str(int(r.countries)), ex_path)
        add(f"{key}_sectors", len(r.sectors.split(",")), str(len(r.sectors.split(","))), ex_path)
        add(f"{key}_code", r.lineage, r.lineage, ex_path)
    bad_pair = vp[(vp.same_lineage_by_alignment) != (vp.levels_shared_v41 >= 5)].iloc[0]
    add("vim_discordant_af", round(float(bad_pair.AF_min), 2), f"{100 * bad_pair.AF_min:.0f}%", vp_path)
    add("vim_discordant_ani", float(bad_pair.ANI_aln), f"{bad_pair.ANI_aln:.3f}%", vp_path)
    add("build_families_median", bs["protein_families_per_plasmid_median"], f"{bs['protein_families_per_plasmid_median']:.0f}", bs_path)
    add("build_kmers_median", bs["kmer_hashes_per_plasmid_median"], str(bs["kmer_hashes_per_plasmid_median"]), bs_path)
    add("build_proteins_predicted", bs["proteins_predicted_all_ids"], f"{bs['proteins_predicted_all_ids']:,}", bs_path)
    add("simpson_L1", dvs["pLIN_L1"], f"{dvs['pLIN_L1']:.4f}", dvs_path)
    add("replicon_types_db", dvs["replicon_types"], str(dvs["replicon_types"]), dvs_path)
    pts = ptu[ptu.set == "confirmatory test set"].set_index("method")
    add("ptu_test_n", int(pts.at["pLIN v4.1 L1", "n"]), str(int(pts.at["pLIN v4.1 L1", "n"])), ptu_path)
    add("ptu_test_ARI_L1", round(pts.at["pLIN v4.1 L1", "ARI"], 2), f"{pts.at['pLIN v4.1 L1', 'ARI']:.2f}", ptu_path)
    add("ptu_test_ARI_MOBprim", round(pts.at["MOB-suite primary", "ARI"], 2), f"{pts.at['MOB-suite primary', 'ARI']:.2f}", ptu_path)

    # tool stability (development data: v4 evaluation plasmids, half coded first, the rest added)
    for (mode, level), key in ((("add", "subcommunity"), "pling_add_sub"), (("add", "community"), "pling_add_comm"),
                               (("rebuild", "subcommunity"), "pling_rebuild_sub"), (("rebuild", "community"), "pling_rebuild_comm")):
        r = ps.loc[(mode, level)]
        add(f"{key}_labels", float(r.label_retained_pct), pct(r.label_retained_pct, 1), ps_path)
        add(f"{key}_splitpct", float(r.pairs_split_pct), pct(r.pairs_split_pct, 1), ps_path)
        add(f"{key}_merged", int(r.pairs_merged), f"{int(r.pairs_merged):,}", ps_path)
    for mode in ("existing", "rebuild"):
        r = ms.loc[mode]
        add(f"mge_{mode}_labels", float(r.label_retained_pct), pct(r.label_retained_pct, 1), ms_path)
        add(f"mge_{mode}_merged", int(r.pairs_merged), f"{int(r.pairs_merged):,}", ms_path)

    # speed and scaling on the same machine (benchmark_speed.py; pLIN v4.1 in speed_v41.json)
    bs2_path = os.path.join(BASE_DIR, "output", "backbone_v4", "speed", "speed_results.tsv")
    bsp = pd.read_csv(bs2_path, sep="\t")
    bsp = bsp[bsp.source == "measured"]
    for tool, key in (("MOB-suite", "mob"), ("pling", "pling")):
        t = bsp[(bsp.tool == tool) & (bsp.step == "type plasmids")].sort_values("n_plasmids")
        for _, r in t.iterrows():
            add(f"speed_{key}_{int(r.n_plasmids)}_s", float(r.wall_s), f"{r.wall_s:,.0f} s", bs2_path)

    # cluster sizes per level (release codes)
    rc_codes = os.path.join(V41, "release", "plin_v41_codes.tsv.gz")
    cc = pd.read_csv(rc_codes, sep="\t", dtype=str)
    for k in range(1, 7):
        sizes = cc.groupby(cc.pLIN_v41.str.split(".").str[:k].str.join(".")).size()
        add(f"singleton_L{k}", round(float((sizes == 1).sum() / len(sizes)), 3), pct(100 * (sizes == 1).sum() / len(sizes)), rc_codes)
        add(f"largest_L{k}", int(sizes.max()), f"{int(sizes.max()):,}", rc_codes)

    # key resistance genes in the largest L6 clones (v41_clone_amr.py)
    ca_path = os.path.join(V41, "onehealth", "clone_amr_matrix.tsv")
    ca = pd.read_csv(ca_path, sep="\t", index_col=0)
    genes = ca.drop(columns="plasmids")
    add("clone_amr_n", len(ca), str(len(ca)), ca_path)
    dom = int((genes.max(axis=1) >= 75).sum())
    add("clone_amr_dominant75", dom, str(dom), ca_path)
    add("clone_amr_genes", genes.shape[1], str(genes.shape[1]), ca_path)

    vc_path = os.path.join(V41, "case_studies", "swiss_vim1_codes.tsv")
    vc = pd.read_csv(vc_path, sep="\t")
    vim = vc[vc.key_genes.fillna("").str.contains("blaVIM-1")]
    add("vim_isolates", vc.plasmid_id.str.split("_").str[0].nunique(), str(vc.plasmid_id.str.split("_").str[0].nunique()), vc_path)
    add("vim_contigs", len(vc), str(len(vc)), vc_path)
    add("vim_size_kb", [int(vim.length_bp.min() // 1000), int(vim.length_bp.max() // 1000)],
        f"{vim.length_bp.min() // 1000} to {vim.length_bp.max() // 1000} kb", vc_path)

    v4e_path = os.path.join(BASE_DIR, "output", "backbone_v4", "evaluation", "primary_endpoints.json")
    v4e = json.load(open(v4e_path))
    e = v4e["endpoint1_related_backbone_F1_v4L3_minus_MOBprimary"]
    add("v4_E1", round(e["difference"], 2), f"+{e['difference']:.2f} (95% CI {e['lo']:.2f} to {e['hi']:.2f})", v4e_path)
    add("v4_test_pairs", v4e["test_pairs"], f"{v4e['test_pairs']:,}", v4e_path)
    cmp_path = os.path.join(V41, "dev", "compare.tsv")
    dv_cmp = pd.read_csv(cmp_path, sep="\t")
    for eng, short in (("founder", "fo"), ("hybrid", "hy"), ("nn", "nn")):
        g = dv_cmp[dv_cmp.engine == eng]
        bb = [float(x.split("@")[0]) for x in g.backbone_F1_best]
        ln = [float(x.split("@")[0]) for x in g.lineage_F1_best]
        add(f"dev_{short}_n", len(g), str(len(g)), cmp_path)
        add(f"dev_{short}_bb", round(max(bb), 2), f"{max(bb):.2f}", cmp_path)
        add(f"dev_{short}_lin", round(max(ln), 2), f"{max(ln):.2f}", cmp_path)
        add(f"dev_{short}_requery", round(g.requery.min(), 2), f"{100 * g.requery.min():.0f}%", cmp_path)

    add("dev_designs", len(pd.read_csv(cmp_path, sep="\t")), str(len(pd.read_csv(cmp_path, sep="\t"))), cmp_path)
    sha_path = os.path.join(V41, "PREREGISTRATION_v4.1.sha256")
    sha = open(sha_path).read().split()[0]
    add("prereg_sha256", sha, sha, sha_path)

    # sensitivity of the confirmatory comparison to the truth definitions (v41_truth_sensitivity.py)
    ts_path = os.path.join(V41, "confirm", "truth_sensitivity_best.tsv")
    ts = pd.read_csv(ts_path, sep="\t")
    best = ts[ts.tool != "registered comparison"].pivot_table(index=["truth", "definition"], columns="tool", values="F1_w")
    others = best.drop(columns=["pLIN v4.1", "pLIN v4"]).max(axis=1)
    lin, bb = best.loc["same_lineage"], best.loc["related_backbone"]
    lin_o, bb_o = others.loc["same_lineage"], others.loc["related_backbone"]
    reg = ts[ts.tool == "registered comparison"]
    lin_reg, bb_reg = reg[reg.truth == "same_lineage"], reg[reg.truth == "related_backbone"]
    add("ts_lineage_defs", len(lin), str(len(lin)), ts_path)
    add("ts_lineage_top", int((lin["pLIN v4.1"] > lin_o).sum()), str(int((lin["pLIN v4.1"] > lin_o).sum())), ts_path)
    add("ts_lineage_plin_range", [round(lin["pLIN v4.1"].min(), 2), round(lin["pLIN v4.1"].max(), 2)],
        f"{lin['pLIN v4.1'].min():.2f} to {lin['pLIN v4.1'].max():.2f}", ts_path)
    add("ts_lineage_other_range", [round(lin_o.min(), 2), round(lin_o.max(), 2)],
        f"{lin_o.min():.2f} to {lin_o.max():.2f}", ts_path)
    add("ts_lineage_reg_sig", int((lin_reg.diff_lo > 0).sum()), str(int((lin_reg.diff_lo > 0).sum())), ts_path)
    add("ts_backbone_defs", len(bb), str(len(bb)), ts_path)
    add("ts_backbone_sig", int((bb_reg.diff_lo > 0).sum()), str(int((bb_reg.diff_lo > 0).sum())), ts_path)
    add("ts_backbone_worse", int((bb_reg.diff_hi < 0).sum()), str(int((bb_reg.diff_hi < 0).sum())), ts_path)
    # external validation on independent hospital datasets (v41_external.py; PREREGISTRATION_external_v4.1.md)
    xe_path = os.path.join(V41, "external", "external_endpoints.json")
    xe = json.load(open(xe_path))
    xm_path = os.path.join(V41, "external", "external_metrics.tsv")
    xm = pd.read_csv(xm_path, sep="\t")
    for ds in ("E1", "E3"):
        d = xe["datasets"][ds]
        add(f"ext_{ds}_n", d["after_excluding_development"], str(d["after_excluding_development"]), xe_path)
        add(f"ext_{ds}_identical", d["in_database_identical"], str(d["in_database_identical"]), xe_path)
    add("ext_n", sum(xe["datasets"][k]["after_excluding_development"] for k in ("E1", "E3")),
        str(sum(xe["datasets"][k]["after_excluding_development"] for k in ("E1", "E3"))), xe_path)
    add("ext_pairs", xe["pairs"], f"{xe['pairs']:,}", xe_path)
    add("ext_lineage_pairs", xe["same_lineage_pairs"], f"{xe['same_lineage_pairs']:,}", xe_path)
    add("ext_backbone_pairs", xe["related_backbone_pairs"], f"{xe['related_backbone_pairs']:,}", xe_path)
    for key, e in (("ext_P1", xe["P1_lineage_L5_vs_MOBsecondary"]), ("ext_P2", xe["P2_backbone_L1_vs_MOBprimary"])):
        add(f"{key}_F1", round(e["F1"], 2), f"{e['F1']:.2f}", xe_path)
        add(f"{key}_F1_CI", [round(e["F1_lo"], 2), round(e["F1_hi"], 2)], f"{e['F1_lo']:.2f} to {e['F1_hi']:.2f}", xe_path)
        add(f"{key}_F1_MOB", round(e["F1_comparator"], 2), f"{e['F1_comparator']:.2f}", xe_path)
        add(f"{key}_diff", round(e["difference"], 2),
            f"{e['difference']:+.2f} (95% CI {e['lo']:.2f} to {e['hi']:.2f})", xe_path)
    f1 = lambda scope, meth, t: xm[(xm.scope == scope) & (xm.method == meth) & (xm.truth == t)].F1_w.iloc[0]
    for scope in ("E1", "E3"):
        for meth, k in (("pLIN v4.1 L5", "L5"), ("MOB-suite secondary", "MOBsec"), ("pLIN v4.1 L1", "L1"),
                        ("MOB-suite primary", "MOBpri")):
            t = "same_lineage" if k in ("L5", "MOBsec") else "related_backbone"
            add(f"ext_{scope}_{k}", round(f1(scope, meth, t), 2), f"{f1(scope, meth, t):.2f}", xm_path)
    for k in ("L2", "L3"):
        add(f"ext_bb_{k}", round(f1("pooled", f"pLIN v4.1 {k}", "related_backbone"), 2),
            f"{f1('pooled', f'pLIN v4.1 {k}', 'related_backbone'):.2f}", xm_path)
    add("ext_pling_lineage", round(f1("pooled", "pling subcommunity", "same_lineage"), 2),
        f"{f1('pooled', 'pling subcommunity', 'same_lineage'):.2f}", xm_path)
    add("ext_pling_backbone", round(f1("pooled", "pling subcommunity", "related_backbone"), 2),
        f"{f1('pooled', 'pling subcommunity', 'related_backbone'):.2f}", xm_path)
    p = xm[(xm.scope == "E3") & (xm.method == "pLIN v4.1 L1") & (xm.truth == "related_backbone")].iloc[0]
    add("ext_E3_L1_precision", round(p.precision_w, 2), f"{p.precision_w:.2f}", xm_path)
    sec = xe["carbapenemase_plasmids_grouped_pct"]
    for ds in ("E1", "E3"):
        add(f"ext_{ds}_carriers", sec[ds]["carrier_plasmids"], str(sec[ds]["carrier_plasmids"]), xe_path)
        add(f"ext_{ds}_carrier_pairs", sec[ds]["carrier_pairs"], f"{sec[ds]['carrier_pairs']:,}", xe_path)
        for meth, k in (("pLIN v4.1 L5", "L5"), ("pLIN v4.1 L6", "L6"), ("MOB-suite secondary", "MOBsec")):
            add(f"ext_{ds}_carb_{k}", sec[ds][meth], f"{sec[ds][meth]:.1f}%", xe_path)
    tb_path = os.path.join(V41, "confirm", "taxon_breakdown.json")
    if os.path.exists(tb_path):
        tb = json.load(open(tb_path))
        for key, short in (("gram_positive", "gp"), ("nonfermenter", "nf"), ("enterobacterales", "ent")):
            if key in tb:
                add(f"taxon_{short}_pairs", tb[key]["pairs"], f"{tb[key]['pairs']:,}", tb_path)
                add(f"taxon_{short}_plin", round(tb[key]["pLIN_L5"], 2), f"{tb[key]['pLIN_L5']:.2f}", tb_path)
                add(f"taxon_{short}_mob", round(tb[key]["MOB_secondary"], 2), f"{tb[key]['MOB_secondary']:.2f}", tb_path)

    mc_path = os.path.join(V41, "confirm", "measure_comparison.json")
    if os.path.exists(mc_path):
        mc = json.load(open(mc_path))
        add("measure_pairs", mc["pairs"], f"{mc['pairs']:,}", mc_path)
        for key, short in (("cosine", "cos"), ("prot_cont", "protc"), ("prot_jacc", "protj"),
                           ("kcont", "kcont"), ("kmin", "kmin")):
            v = mc["measures"][key]
            add(f"meas_{short}_lin", v["AUC_same_lineage"], f"{v['AUC_same_lineage']:.2f}", mc_path)
            add(f"meas_{short}_bb", v["AUC_related_backbone"], f"{v['AUC_related_backbone']:.2f}", mc_path)

    cs_path = os.path.join(V41, "confirm", "comparator_sensitivity.tsv")
    cs = pd.read_csv(cs_path, sep="\t")
    n_settings = cs.assign(t=cs.tool.str.split().str[0]).drop_duplicates(["t", "setting"]).shape[0]   # pling and mge-cluster runs
    add("cs_settings", n_settings, str(n_settings), cs_path)
    for truth, key in (("same_lineage", "lineage"), ("related_backbone", "backbone")):
        d = cs[cs.truth == truth]
        b = d.loc[d.F1_w.idxmax()]
        add(f"cs_{key}_best", round(b.F1_w, 2), f"{b.F1_w:.2f}", cs_path)
        add(f"cs_{key}_best_setting", f"{b.tool}, {b.setting}", f"{b.tool} with {b.setting.replace('dcj', 'DCJ-Indel threshold')}", cs_path)
        mg = d[d.tool == "mge-cluster"].F1_w.max()
        add(f"cs_{key}_mge_max", round(mg, 2), f"{mg:.2f}", cs_path)
    for t in (0.3, 0.4, 0.7):
        k = f"backbone AF>={t}"
        add(f"ts_bb{int(t * 10)}_plin", round(bb.loc[k, "pLIN v4.1"], 2), f"{bb.loc[k, 'pLIN v4.1']:.2f}", ts_path)
        add(f"ts_bb{int(t * 10)}_other", round(bb_o.loc[k], 2), f"{bb_o.loc[k]:.2f}", ts_path)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(facts, open(OUT, "w"), indent=1)
    for k, v in facts.items():
        print(f"{k:32s} {v['text']}")


if __name__ == "__main__":
    main()
