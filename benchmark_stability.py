#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Label stability as the database grows: pLIN vs pling (vs MOB-suite).

A typing scheme is only useful as a nomenclature if a plasmid keeps its label
when new plasmids are added. We simulate database growth on the 2,747
comparator-benchmark plasmids: snapshot A is a random subset, the grown
database is all 2,747 (A + B). For the plasmids in A we compare the label they
had in snapshot A with the label they have after growth.

Metrics (over plasmids in A; unassigned plasmids count as singletons)
  label_retained_pct   share of A plasmids whose label string is unchanged
  pairs_split_pct      of A pairs grouped together before, share now separated
  pairs_merged         A pairs separated before that are now grouped together
  ARI                  adjusted Rand index between before and after partitions

Modes
  pLIN incremental   tree built on A, then B added in accession order; A plasmids
                     re-queried against the grown tree (intended use)
  pLIN rebuild       codes rebuilt from scratch on A + B (not the intended use)
  pling rebuild      pling cluster on A, and pling cluster on A + B
  pling add          pling cluster on A, then pling add B (pling's update mode)
  MOB-suite          mob_typer types each plasmid against its fixed reference
                     database, so labels do not depend on the other query
                     plasmids; they change only when the reference database is
                     re-released. Not simulated here.

pLIN is cheap, so it is run on 5 random seeds and snapshot sizes 25/50/75%;
pling (hours per run) on seed 1 at 50% only.

Usage:
  python benchmark_stability.py
Output: output/comparator_benchmark/stability/{stability_results.tsv, plin_codes_seed1.tsv}
"""

import os
import random
import warnings

import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.metrics import adjusted_rand_score

from plin_founder import PLIN_THRESHOLDS, FounderTree, build_founder_codes, kmer_vector

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CB = os.path.join(BASE_DIR, "output", "comparator_benchmark")
ST = os.path.join(CB, "stability")
warnings.filterwarnings("ignore", message="The number of unique classes")  # sklearn, harmless for ARI
SEEDS = [1, 2, 3, 4, 5]
FRACTIONS = [0.25, 0.5, 0.75]


def plasmid_id(path):
    return os.path.basename(path)[:-len(".fasta")]


def snapshot(paths, frac, seed):
    """Same sampling as the pling snapshot files (seed 1, 50% -> snapshotA_seed1.txt)."""
    return sorted(random.Random(seed).sample(paths, int(len(paths) * frac)))


def compare(before, after):
    """Stability metrics for two labelings (dicts id -> label, None = unassigned) over the same ids."""
    ids = sorted(before)
    b = [before[i] if before[i] is not None else f"_single_{i}" for i in ids]
    a = [after.get(i) if after.get(i) is not None else f"_single_{i}" for i in ids]
    tab = pd.crosstab(pd.Series(b, name="b"), pd.Series(a, name="a")).values
    c2 = lambda x: x * (x - 1) // 2
    both = int(c2(tab).sum())
    tog_before, tog_after = int(c2(tab.sum(1)).sum()), int(c2(tab.sum(0)).sum())
    return {"n_plasmids": len(ids),
            "label_retained_pct": round(100 * np.mean([x == y for x, y in zip(b, a)]), 2),
            "pairs_together_before": tog_before,
            "pairs_split_pct": round(100 * (tog_before - both) / tog_before, 2) if tog_before else np.nan,
            "pairs_merged": tog_after - both,
            "ARI": round(adjusted_rand_score(b, a), 4),
            "unassigned_before": sum(before[i] is None for i in ids),
            "unassigned_after": sum(after.get(i) is None for i in ids)}


def prefixes(codes, k):
    return {i: ".".join(map(str, c[:k])) for i, c in codes.items()}


def plin_runs(paths, vec):
    rows, saved = [], None
    for frac in FRACTIONS:
        for seed in SEEDS:
            A = [plasmid_id(p) for p in snapshot(paths, frac, seed)]
            setA = set(A)
            B = sorted(i for i in vec if i not in setA)
            codesA_list, tree = build_founder_codes([vec[i] for i in A], A, PLIN_THRESHOLDS)
            codesA = dict(zip(A, codesA_list))
            for i in B:                                   # incremental growth
                tree.assign(vec[i], key=i)
            requery = {i: tree.assign(vec[i], key=i, add=False)[0] for i in A}
            allids = sorted(vec)
            rebuilt = dict(zip(allids, build_founder_codes([vec[i] for i in allids], allids, PLIN_THRESHOLDS)[0]))
            for mode, after in (("incremental", requery), ("rebuild", rebuilt)):
                for k in range(1, 7):
                    rows.append({"tool": "pLIN", "mode": mode, "level": f"L{k}", "snapshot_frac": frac, "seed": seed,
                                 **compare(prefixes(codesA, k), prefixes(after, k))})
            if frac == 0.5 and seed == 1:
                saved = pd.DataFrame({"plasmid_id": A,
                                      "pLIN_snapshotA": [".".join(map(str, codesA[i])) for i in A],
                                      "pLIN_incremental": [".".join(map(str, requery[i])) for i in A],
                                      "pLIN_rebuild": [".".join(map(str, rebuilt[i])) for i in A]})
    return rows, saved


def load_pling(out_dir):
    """{'community': {...}, 'subcommunity': {...}} or None if the run has not finished."""
    import glob
    typing = glob.glob(os.path.join(out_dir, "dcj_thresh_*_graph", "objects", "typing.tsv"))
    if not typing:
        return None
    files = {"community": os.path.join(out_dir, "containment", "containment_communities", "objects", "communities.tsv"),
             "subcommunity": typing[0]}
    return {col: dict(pd.read_csv(f, sep="\t", dtype=str).iloc[:, :2].values) for col, f in files.items()}


def pling_runs(A):
    rows = []
    before = load_pling(os.path.join(ST, "pling_A"))
    if before is None:
        print("pling snapshot A not finished; skipping pling")
        return rows
    for mode, d in (("rebuild", os.path.join(CB, "pling_out")), ("add", os.path.join(ST, "pling_A_add"))):
        after = load_pling(d)
        if after is None:
            print(f"pling {mode} run not finished; skipping")
            continue
        for level in ("community", "subcommunity"):
            rows.append({"tool": "pling", "mode": mode, "level": level, "snapshot_frac": 0.5, "seed": 1,
                         **compare({i: before[level].get(i) for i in A}, after[level])})
    return rows


def main():
    paths = [l.strip() for l in open(os.path.join(CB, "pling_inputs.txt")) if l.strip()]
    cache = os.path.join(ST, "vectors.npz")
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        vec = dict(zip(z["ids"].tolist(), z["X"]))
    else:
        vec = {plasmid_id(p): kmer_vector(str(next(SeqIO.parse(p, "fasta")).seq)) for p in paths}
        np.savez_compressed(cache, ids=np.array(list(vec)), X=np.array(list(vec.values())))

    rows, saved = plin_runs(paths, vec)
    saved.to_csv(os.path.join(ST, "plin_codes_seed1.tsv"), sep="\t", index=False)
    A = [plasmid_id(l.strip()) for l in open(os.path.join(ST, "snapshotA_seed1.txt")) if l.strip()]
    assert A == [plasmid_id(p) for p in snapshot(paths, 0.5, 1)], "pling snapshot file does not match seed-1 sampling"
    rows += pling_runs(A)

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(ST, "stability_results.tsv"), sep="\t", index=False)
    metrics = ["label_retained_pct", "pairs_split_pct", "pairs_merged", "ARI"]
    summ = res.groupby(["tool", "mode", "level", "snapshot_frac"])[metrics].agg(["mean", "min"]).round(3)
    pd.set_option("display.width", 250)
    print(summ.to_string())


if __name__ == "__main__":
    main()
