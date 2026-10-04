#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Alignment-based validation of pLIN codes beyond ANI: alignment fraction,
backbone vs accessory coverage, and the actual sequence differences between
plasmids that share an identical or near-identical pLIN code.

Addresses reviewer requests (JOI R1 major 1; npj AMR R2) that validation go
beyond FastANI, which reports identity only over the fragments it can map and
so can stay high when two modular plasmids share little of their architecture.

Design
  1. Pairs are stratified by the deepest pLIN level they share
     (nested single-linkage, so "deepest shared level" is well defined):
        L6  identical six-level code
        L5  share L1-L5, differ at L6   ("near-identical" code)
        L4  share L1-L4, differ at L5
        L3  share L1-L3, differ at L4   (coarse; reference stratum)
     Within each stratum, up to --pairs-per-group pairs are drawn per replicon
     group (inc_type of the first member). Clusters are sampled uniformly
     before members, so a few very large clusters do not dominate.
  2. Each pair is aligned with blastn (megablast, e-value <= 1e-10, hits
     >= 100 bp). Hits are reduced to a one-to-one set by greedy selection on
     bitscore, rejecting any hit that overlaps an accepted hit by >10% of its
     length on either sequence (repeat copies, e.g. multiple IS26, are
     otherwise counted many times).
  3. From the one-to-one set, per pair:
        AF_shorter / AF_longer  fraction of each plasmid covered
        AF_min                  min of the two (bidirectional AF)
        ANI_aln                 length-weighted identity over aligned blocks
        SNPs_per_100kb          mismatches per 100 kb aligned
        n_blocks_1kb            aligned blocks >= 1 kb (fragmentation)
  4. Each plasmid is partitioned into accessory (AMRFinderPlus AMR elements,
     curated IS-element hits, curated composite-transposon spans) and
     backbone (everything else). AF is recomputed separately for backbone
     and accessory bases, and the unaligned sequence of each plasmid is
     decomposed into backbone vs accessory bases. This tests directly
     whether shared codes reflect a shared backbone or shared cargo.
  5. Worked examples: for the L6 and L5 strata, the pair closest to the
     stratum median AF_min is chosen in each of three taxonomic categories
     (Enterobacterales, Gram-positive, non-fermenter), plus the lowest-AF_min
     L6 pair as the worst case. Selection is rule-based, not hand-picked.
  6. Chaining diagnostics. Single linkage joins two plasmids at a level if
     any chain of members links them, so a shared prefix need not mean the
     two are directly within that level's threshold. Each pair is flagged
     for (i) membership of the largest cluster at its shared level and
     (ii) whether its direct cosine distance is within that threshold.

Usage:
  python validate_alignment_backbone.py
  python validate_alignment_backbone.py --pairs-per-group 15 --threads 12 --seed 42
  python validate_alignment_backbone.py --assignments output/scheme_benchmark/assignments_founder_global.tsv \
      --out-dir output/scheme_benchmark/alignment_founder_global

Output (output/alignment_validation/):
  pair_metrics.tsv        one row per pair
  stratum_summary.tsv     median [IQR] per stratum
  stats.json              tests reported in the manuscript
  worked_examples.tsv     selected example pairs + their differences
  example_blocks/*.tsv    one-to-one blocks for each worked example
"""

import os
import json
import glob
import random
import argparse
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import spearmanr, kruskal, wilcoxon

from plin_founder import PLIN_THRESHOLDS, kmer_vector

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TRAINING_DIR = os.path.join(BASE_DIR, "plasmid_sequences_for_training")
ASSIGN_TSV = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")
AMR_TSV = os.path.join(BASE_DIR, "output", "amrfinder", "amrfinder_all_plasmids.tsv")
IS_TSV = os.path.join(BASE_DIR, "output", "mge_detection", "is_element_hits_curated.tsv")
TN_TSV = os.path.join(BASE_DIR, "output", "mge_detection", "composite_transposons_curated.tsv")
OUT_DIR = os.path.join(BASE_DIR, "output", "alignment_validation")

BINS = ["bin_A", "bin_B", "bin_C", "bin_D", "bin_E", "bin_F"]
STRATA = [6, 5, 4, 3]
# cosine-distance threshold of each level
LEVEL_THRESHOLD = {i + 1: t for i, t in enumerate(PLIN_THRESHOLDS.values())}
MIN_HIT_LEN = 100
MAX_OVERLAP = 0.10

GRAM_POSITIVE = {"repSA_large", "repSA_small", "repEF_conj", "repEF_res"}
NON_FERMENTER = {"repAci1", "repAci_large", "repPae_large", "repPae_small"}


def category(inc):
    if inc in GRAM_POSITIVE:
        return "Gram-positive"
    if inc in NON_FERMENTER:
        return "Non-fermenter"
    return "Enterobacterales"


# ── Inputs ────────────────────────────────────────────────────────────────────

def fasta_paths():
    """plasmid_id -> fasta path (first occurrence; 21 E. faecium IDs are in two folders)."""
    paths = {}
    for p in sorted(glob.glob(os.path.join(TRAINING_DIR, "*", "fastas", "*.fasta"))):
        paths.setdefault(os.path.basename(p)[:-len(".fasta")], p)
    return paths


def load_assignments():
    df = pd.read_csv(ASSIGN_TSV, sep="\t")
    return df.drop_duplicates("plasmid_id").reset_index(drop=True)


def load_accessory_intervals():
    """plasmid_id -> list of (start, end) 1-based closed accessory intervals, plus AMR gene sets."""
    intervals, amr_genes = {}, {}

    amr = pd.read_csv(AMR_TSV, sep="\t", low_memory=False)
    amr = amr[amr["Type"] == "AMR"]
    for pid, start, stop, sym in zip(amr["Contig id"], amr["Start"], amr["Stop"], amr["Element symbol"]):
        intervals.setdefault(pid, []).append((int(start), int(stop)))
        amr_genes.setdefault(pid, set()).add(sym)

    is_hits = pd.read_csv(IS_TSV, sep="\t")
    for pid, s, e in zip(is_hits["source_file"], is_hits["q_start"], is_hits["q_end"]):
        intervals.setdefault(pid, []).append((min(s, e), max(s, e)))

    tn = pd.read_csv(TN_TSV, sep="\t")
    for pid, a, b in zip(tn["source_file"], tn["is1_pos"], tn["is2_pos"]):
        coords = [int(x) for x in f"{a}-{b}".split("-")]
        intervals.setdefault(pid, []).append((min(coords), max(coords)))

    return intervals, amr_genes


def cosine_distance(u, v):
    return 1.0 - float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))


# ── Pair sampling ─────────────────────────────────────────────────────────────

def sample_pairs(df, pairs_per_group, seed):
    """Sample pairs stratified by deepest shared pLIN level and by replicon group."""
    rng = random.Random(seed)
    keys = df[BINS].astype(str).agg(".".join, axis=1)
    prefixes = {lvl: df[BINS[:lvl]].astype(str).agg(".".join, axis=1) for lvl in range(1, 7)}
    largest = {lvl: prefixes[lvl].value_counts().idxmax() for lvl in range(1, 7)}
    pairs, seen = [], set()

    for lvl in STRATA:
        pref = prefixes[lvl]
        # child level: pair members must differ here (none for L6)
        child = prefixes[lvl + 1] if lvl < 6 else None
        for inc, sub in df.groupby("inc_type"):
            clusters = [c for c, g in pref[sub.index].groupby(pref[sub.index]) if True]
            # a cluster is usable if it contains >=2 plasmids (any inc type) that
            # differ at the child level (or are identical-code when lvl == 6)
            usable = []
            for c in clusters:
                members = df.index[pref == c].tolist()
                anchors = [i for i in members if i in sub.index]
                if lvl == 6:
                    ok = len(members) >= 2
                else:
                    ok = len(set(child[members])) >= 2
                if ok and anchors:
                    usable.append((c, members, anchors))
            if not usable:
                continue
            attempts, got = 0, 0
            while got < pairs_per_group and attempts < pairs_per_group * 50:
                attempts += 1
                c, members, anchors = rng.choice(usable)
                a = rng.choice(anchors)
                partners = [m for m in members if m != a and
                            (lvl == 6 or child[m] != child[a])]
                if not partners:
                    continue
                b = rng.choice(partners)
                key = (lvl, min(a, b), max(a, b))
                if key in seen:
                    continue
                seen.add(key)
                pairs.append({"stratum": f"L{lvl}", "shared_levels": lvl, "a": a, "b": b,
                              "anchor_inc": inc, "anchor_cluster": c,
                              "cluster_size": len(members),
                              "in_largest_cluster": c == largest[lvl]})
                got += 1
    return pairs


# ── Alignment ─────────────────────────────────────────────────────────────────

def run_blast(path_q, path_s):
    cmd = ["blastn", "-query", path_q, "-subject", path_s, "-evalue", "1e-10",
           "-outfmt", "6 qstart qend sstart send pident length mismatch gapopen bitscore",
           "-max_hsps", "1000", "-dust", "no"]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    hits = []
    for line in out.splitlines():
        qs, qe, ss, se, pid, ln, mm, go, bs = line.split("\t")
        ln = int(ln)
        if ln < MIN_HIT_LEN:
            continue
        hits.append({"qs": int(qs), "qe": int(qe),
                     "ss": min(int(ss), int(se)), "se": max(int(ss), int(se)),
                     "strand": "+" if int(se) >= int(ss) else "-",
                     "pident": float(pid), "length": ln, "mismatch": int(mm),
                     "gapopen": int(go), "bitscore": float(bs)})
    return hits


def overlap_len(s, e, accepted):
    total = 0
    for a, b in accepted:
        lo, hi = max(s, a), min(e, b)
        if hi >= lo:
            total += hi - lo + 1
    return total


def one_to_one(hits):
    """Greedy one-to-one hit set: best bitscore first, reject >10% overlap on either sequence."""
    acc, q_iv, s_iv = [], [], []
    for h in sorted(hits, key=lambda h: -h["bitscore"]):
        qlen = h["qe"] - h["qs"] + 1
        slen = h["se"] - h["ss"] + 1
        if overlap_len(h["qs"], h["qe"], q_iv) > MAX_OVERLAP * qlen:
            continue
        if overlap_len(h["ss"], h["se"], s_iv) > MAX_OVERLAP * slen:
            continue
        acc.append(h)
        q_iv.append((h["qs"], h["qe"]))
        s_iv.append((h["ss"], h["se"]))
    return acc


def mask_from_intervals(length, intervals):
    m = np.zeros(length + 1, dtype=bool)   # 1-based; index 0 unused
    for s, e in intervals:
        m[max(1, s):min(length, e) + 1] = True
    return m


def align_pair(task):
    """Worker: align one pair and compute AF / backbone metrics."""
    pa, pb, la, lb, acc_a, acc_b, keep_blocks = task
    hits = one_to_one(run_blast(pa, pb))

    cov_a = mask_from_intervals(la, [(h["qs"], h["qe"]) for h in hits])
    cov_b = mask_from_intervals(lb, [(h["ss"], h["se"]) for h in hits])
    accm_a = mask_from_intervals(la, acc_a)
    accm_b = mask_from_intervals(lb, acc_b)

    def split(cov, accm, length):
        cov, accm = cov[1:], accm[1:]
        bb, ac = ~accm, accm
        return {
            "AF": cov.sum() / length,
            "bb_len": int(bb.sum()), "acc_len": int(ac.sum()),
            "bb_AF": (cov & bb).sum() / bb.sum() if bb.sum() else np.nan,
            "acc_AF": (cov & ac).sum() / ac.sum() if ac.sum() else np.nan,
            "unaligned": int((~cov).sum()),
            "unaligned_acc": int((~cov & ac).sum()),
        }

    ra, rb = split(cov_a, accm_a, la), split(cov_b, accm_b, lb)
    aln_len = sum(h["length"] for h in hits)
    ani = sum(h["pident"] * h["length"] for h in hits) / aln_len if aln_len else np.nan

    # identity over backbone-only portion of each block (block pident weighted by its backbone overlap on A)
    bb_w = [(h["pident"], int((~accm_a[h["qs"]:h["qe"] + 1]).sum())) for h in hits]
    bb_tot = sum(w for _, w in bb_w)
    ani_bb = sum(p * w for p, w in bb_w) / bb_tot if bb_tot else np.nan

    res = {"a": ra, "b": rb, "n_hits": len(hits),
           "n_blocks_1kb": sum(h["length"] >= 1000 for h in hits),
           "largest_block": max((h["length"] for h in hits), default=0),
           "aligned_bp": aln_len, "ANI_aln": ani, "ANI_backbone": ani_bb,
           "mismatches": sum(h["mismatch"] for h in hits),
           "gap_opens": sum(h["gapopen"] for h in hits),
           "n_inverted_blocks": sum(h["strand"] == "-" and h["length"] >= 1000 for h in hits)}
    if keep_blocks:
        res["blocks"] = hits
    return res


# ── Main ──────────────────────────────────────────────────────────────────────

def iqr_str(x):
    x = pd.Series(x).dropna()
    if x.empty:
        return "NA"
    return f"{x.median():.3f} [{x.quantile(.25):.3f}–{x.quantile(.75):.3f}]"


def main():
    global ASSIGN_TSV, OUT_DIR
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs-per-group", type=int, default=15)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--assignments", default=ASSIGN_TSV,
                    help="pLIN assignment table to validate (default: output/pLIN_assignments.tsv)")
    ap.add_argument("--out-dir", default=OUT_DIR)
    args = ap.parse_args()
    ASSIGN_TSV, OUT_DIR = args.assignments, args.out_dir
    os.makedirs(os.path.join(OUT_DIR, "example_blocks"), exist_ok=True)

    print("[1/5] Loading assignments, sequences, annotations ...")
    df = load_assignments()
    paths = fasta_paths()
    intervals, amr_genes = load_accessory_intervals()

    print("[2/5] Sampling pairs ...")
    pairs = sample_pairs(df, args.pairs_per_group, args.seed)
    print(f"  {len(pairs)} pairs: " +
          ", ".join(f"{s}={sum(p['stratum'] == s for p in pairs)}" for s in [f'L{l}' for l in STRATA]))

    used = sorted({p["a"] for p in pairs} | {p["b"] for p in pairs})
    seqs, vecs = {}, {}
    for i in used:
        pid = df.at[i, "plasmid_id"]
        rec = next(SeqIO.parse(paths[pid], "fasta"))
        seqs[i] = str(rec.seq)
        vecs[i] = kmer_vector(seqs[i])

    print(f"[3/5] Aligning {len(pairs)} pairs with blastn on {args.threads} workers ...")
    tasks = []
    for p in pairs:
        a, b = p["a"], p["b"]
        # query = shorter plasmid, so AF_shorter is always side A
        if len(seqs[a]) > len(seqs[b]):
            a, b = b, a
            p["a"], p["b"] = a, b
        pa, pb = df.at[a, "plasmid_id"], df.at[b, "plasmid_id"]
        tasks.append((paths[pa], paths[pb], len(seqs[a]), len(seqs[b]),
                      intervals.get(pa, []), intervals.get(pb, []), p["stratum"] in ("L6", "L5")))
    with ProcessPoolExecutor(max_workers=args.threads) as ex:
        results = list(ex.map(align_pair, tasks, chunksize=4))

    print("[4/5] Assembling metrics ...")
    rows, blocks_by_row = [], {}
    for idx, (p, r) in enumerate(zip(pairs, results)):
        a, b = p["a"], p["b"]
        pa, pb = df.at[a, "plasmid_id"], df.at[b, "plasmid_id"]
        ga, gb = amr_genes.get(pa, set()), amr_genes.get(pb, set())
        union = ga | gb
        rows.append({
            "stratum": p["stratum"], "shared_levels": p["shared_levels"],
            "plasmid_shorter": pa, "plasmid_longer": pb,
            "inc_shorter": df.at[a, "inc_type"], "inc_longer": df.at[b, "inc_type"],
            "same_inc": df.at[a, "inc_type"] == df.at[b, "inc_type"],
            "category": category(p["anchor_inc"]), "anchor_inc": p["anchor_inc"],
            "pLIN_shorter": df.at[a, "pLIN"], "pLIN_longer": df.at[b, "pLIN"],
            "len_shorter": len(seqs[a]), "len_longer": len(seqs[b]),
            "length_ratio": len(seqs[b]) / len(seqs[a]),
            "cosine_distance": cosine_distance(vecs[a], vecs[b]),
            "cluster_size": p["cluster_size"], "in_largest_cluster": p["in_largest_cluster"],
            "direct_within_threshold": cosine_distance(vecs[a], vecs[b]) <= LEVEL_THRESHOLD[p["shared_levels"]],
            "AF_shorter": r["a"]["AF"], "AF_longer": r["b"]["AF"],
            "AF_min": min(r["a"]["AF"], r["b"]["AF"]),
            "ANI_aln": r["ANI_aln"], "ANI_backbone": r["ANI_backbone"],
            "aligned_bp": r["aligned_bp"],
            "SNPs_per_100kb": 1e5 * r["mismatches"] / r["aligned_bp"] if r["aligned_bp"] else np.nan,
            "gap_opens": r["gap_opens"], "n_blocks_1kb": r["n_blocks_1kb"],
            "largest_block": r["largest_block"], "n_inverted_blocks": r["n_inverted_blocks"],
            "bb_AF_shorter": r["a"]["bb_AF"], "acc_AF_shorter": r["a"]["acc_AF"],
            "bb_AF_longer": r["b"]["bb_AF"], "acc_AF_longer": r["b"]["acc_AF"],
            "acc_frac_shorter": r["a"]["acc_len"] / len(seqs[a]),
            "acc_frac_longer": r["b"]["acc_len"] / len(seqs[b]),
            "unaligned_shorter": r["a"]["unaligned"], "unaligned_acc_shorter": r["a"]["unaligned_acc"],
            "unaligned_longer": r["b"]["unaligned"], "unaligned_acc_longer": r["b"]["unaligned_acc"],
            "amr_shorter": ";".join(sorted(ga)), "amr_longer": ";".join(sorted(gb)),
            "amr_jaccard": len(ga & gb) / len(union) if union else np.nan,
            "amr_only_shorter": ";".join(sorted(ga - gb)), "amr_only_longer": ";".join(sorted(gb - ga)),
        })
        if "blocks" in r:
            blocks_by_row[idx] = r["blocks"]
    m = pd.DataFrame(rows)
    # share of the longer plasmid's unaligned sequence that is annotated cargo
    m["unaligned_acc_share_longer"] = np.where(
        m["unaligned_longer"] > 0, m["unaligned_acc_longer"] / m["unaligned_longer"].replace(0, np.nan), np.nan)
    m.to_csv(os.path.join(OUT_DIR, "pair_metrics.tsv"), sep="\t", index=False)

    print("[5/5] Summaries, tests, worked examples ...")
    summ = []
    for s in [f"L{l}" for l in STRATA]:
        g = m[m.stratum == s]
        summ.append({
            "stratum": s, "n_pairs": len(g), "n_replicon_groups": g.anchor_inc.nunique(),
            "pct_same_inc": 100 * g.same_inc.mean(),
            "cosine_distance": iqr_str(g.cosine_distance),
            "AF_shorter": iqr_str(g.AF_shorter), "AF_longer": iqr_str(g.AF_longer),
            "AF_min": iqr_str(g.AF_min), "ANI_aln": iqr_str(g.ANI_aln),
            "SNPs_per_100kb": iqr_str(g.SNPs_per_100kb),
            "bb_AF_shorter": iqr_str(g.bb_AF_shorter), "acc_AF_shorter": iqr_str(g.acc_AF_shorter),
            "length_ratio": iqr_str(g.length_ratio),
            "pct_AFmin_ge_0.8": 100 * (g.AF_min >= 0.8).mean(),
            "pct_AFshorter_ge_0.9_ANI_ge_99": 100 * ((g.AF_shorter >= 0.9) & (g.ANI_aln >= 99)).mean(),
            "pct_AFmin_lt_0.5": 100 * (g.AF_min < 0.5).mean(),
            "amr_jaccard": iqr_str(g.amr_jaccard),
        })
    summ = pd.DataFrame(summ)
    chain = []
    for s_ in [f"L{l}" for l in STRATA]:
        for flag_col in ["in_largest_cluster", "direct_within_threshold"]:
            for val in [True, False]:
                g = m[(m.stratum == s_) & (m[flag_col] == val)]
                chain.append({"stratum": s_, "split": f"{flag_col}={val}", "n_pairs": len(g),
                              "AF_min": iqr_str(g.AF_min), "AF_shorter": iqr_str(g.AF_shorter),
                              "ANI_aln": iqr_str(g.ANI_aln),
                              "pct_AFmin_ge_0.8": 100 * (g.AF_min >= 0.8).mean() if len(g) else np.nan})
    chain = pd.DataFrame(chain)
    chain.to_csv(os.path.join(OUT_DIR, "chaining_summary.tsv"), sep="\t", index=False)
    print(chain.to_string(index=False))
    sizes = {}
    for lvl in range(1, 7):
        vc = df[BINS[:lvl]].astype(str).agg(".".join, axis=1).value_counts()
        sizes[f"L{lvl}"] = {"n_clusters": int(len(vc)), "largest": int(vc.iloc[0]),
                            "pct_in_largest": float(100 * vc.iloc[0] / len(df)),
                            "n_singletons": int((vc == 1).sum())}
    summ.to_csv(os.path.join(OUT_DIR, "stratum_summary.tsv"), sep="\t", index=False)
    print(summ.to_string(index=False))

    stats = {"n_pairs": int(len(m)), "seed": args.seed, "pairs_per_group": args.pairs_per_group,
             "cluster_sizes": sizes}
    groups = [m[m.stratum == s] for s in [f"L{l}" for l in STRATA]]
    for col in ["AF_min", "AF_shorter", "ANI_aln"]:
        h, p = kruskal(*[g[col].dropna() for g in groups])
        stats[f"kruskal_{col}"] = {"H": float(h), "p": float(p)}
    rho, p = spearmanr(m.shared_levels, m.AF_min)
    stats["spearman_sharedlevels_AFmin"] = {"rho": float(rho), "p": float(p)}
    rho, p = spearmanr(m.cosine_distance, m.AF_min)
    stats["spearman_cosine_AFmin"] = {"rho": float(rho), "p": float(p)}
    rho, p = spearmanr(m.cosine_distance, m.ANI_aln, nan_policy="omit")
    stats["spearman_cosine_ANI"] = {"rho": float(rho), "p": float(p)}
    # backbone vs accessory coverage of the shorter plasmid, paired, per stratum
    for s in [f"L{l}" for l in STRATA]:
        g = m[(m.stratum == s)].dropna(subset=["bb_AF_shorter", "acc_AF_shorter"])
        if len(g) >= 10:
            w, p = wilcoxon(g.bb_AF_shorter, g.acc_AF_shorter)
            stats[f"wilcoxon_bb_vs_acc_{s}"] = {
                "n": int(len(g)), "median_bb": float(g.bb_AF_shorter.median()),
                "median_acc": float(g.acc_AF_shorter.median()), "W": float(w), "p": float(p),
                "pct_bb_gt_acc": float(100 * (g.bb_AF_shorter > g.acc_AF_shorter).mean())}
    # what fills the unaligned sequence of the longer plasmid, vs its overall cargo share
    for s in ["L6", "L5"]:
        g = m[(m.stratum == s) & (m.unaligned_longer >= 1000)]
        stats[f"unaligned_composition_{s}"] = {
            "n_pairs_with_>=1kb_unaligned": int(len(g)),
            "median_acc_share_of_unaligned": float(g.unaligned_acc_share_longer.median()) if len(g) else None,
            "median_acc_share_of_whole_plasmid": float(g.acc_frac_longer.median()) if len(g) else None}
    with open(os.path.join(OUT_DIR, "stats.json"), "w") as fh:
        json.dump(stats, fh, indent=2)
    print(json.dumps(stats, indent=2))

    # worked examples: rule-based selection
    ex_idx = []
    for s in ["L6", "L5"]:
        g = m[m.stratum == s]
        med = g.AF_min.median()
        for cat in ["Enterobacterales", "Gram-positive", "Non-fermenter"]:
            gc = g[g.category == cat]
            if len(gc):
                ex_idx.append(((gc.AF_min - med).abs()).idxmin())
    worst = m[m.stratum == "L6"].AF_min.idxmin()
    if worst not in ex_idx:
        ex_idx.append(worst)
    ex = m.loc[ex_idx].copy()
    ex["example_rule"] = ["closest to stratum median AF_min"] * (len(ex) - 1) + ["lowest AF_min in L6 (worst case)"]
    ex.to_csv(os.path.join(OUT_DIR, "worked_examples.tsv"), sep="\t", index=True, index_label="pair_row")
    for i in ex_idx:
        pd.DataFrame(blocks_by_row[i]).to_csv(
            os.path.join(OUT_DIR, "example_blocks", f"pair{i}_{m.at[i, 'plasmid_shorter']}__{m.at[i, 'plasmid_longer']}.tsv"),
            sep="\t", index=False)
    print(ex[["stratum", "category", "plasmid_shorter", "plasmid_longer", "AF_shorter", "AF_longer",
              "ANI_aln", "SNPs_per_100kb", "amr_only_longer"]].to_string())
    print(f"\nDone. Outputs in {OUT_DIR}")


if __name__ == "__main__":
    main()
