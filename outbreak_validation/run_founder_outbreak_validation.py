#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Outbreak validation of founder-based pLIN codes (74 plasmids, 26 studies).

Replaces, for founder-assigned codes, the nearest-neighbour leave-one-out
benchmark (benchmark_outbreak_detection.py), which applied to the earlier
single-linkage codes.

For each outbreak plasmid (outbreak_validation/*.fasta and
outbreak_validation/expanded_sequences/*.fasta, metadata from
output/outbreak_validation_combined_results.tsv):

  1. Release code: query against the release founder tree
     (data/plin_founder_tree_reference.npz) with the same rule that built it.
  2. Leave-one-out (LOO), for plasmids that are themselves in the database:
     rebuild the database without the plasmid, query it, and test whether its
     cluster-mates at each level are the same as in the release. A plasmid
     that founded no cluster cannot change the tree, so its LOO code equals its
     release code by construction; only founders need a rebuild.
  3. Outbreak grouping: for every pair of outbreak plasmids, whether they
     share an L6 / L5 / L3 code, split by same study vs different study.

Usage:
  python outbreak_validation/run_founder_outbreak_validation.py
Output:
  output/outbreak_validation_founder_results.tsv
  output/outbreak_validation_founder_summary.json
"""

import os
import sys
import glob
import json
from itertools import combinations

import numpy as np
import pandas as pd
from Bio import SeqIO

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BASE_DIR)
from plin_founder import FounderTree, PLIN_THRESHOLDS, kmer_vector  # noqa: E402

META_TSV = os.path.join(BASE_DIR, "output", "outbreak_validation_combined_results.tsv")
REF_TSV = os.path.join(BASE_DIR, "output", "pLIN_reference_assignments.tsv")
TRAIN_TSV = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")
TRAIN_TREE = os.path.join(BASE_DIR, "data", "plin_founder_tree_training.npz")
REF_TREE = os.path.join(BASE_DIR, "data", "plin_founder_tree_reference.npz")
REF_VECTORS = os.path.join(BASE_DIR, "output", "reference_kmer_vectors.npz")
TRAINING_DIR = os.path.join(BASE_DIR, "plasmid_sequences_for_training")
OUT_TSV = os.path.join(BASE_DIR, "output", "outbreak_validation_founder_results.tsv")
OUT_JSON = os.path.join(BASE_DIR, "output", "outbreak_validation_founder_summary.json")
BINS = [f"bin_{l}" for l in PLIN_THRESHOLDS]


def outbreak_fasta(acc):
    for d in (os.path.join(BASE_DIR, "outbreak_validation"),
              os.path.join(BASE_DIR, "outbreak_validation", "expanded_sequences")):
        p = os.path.join(d, f"{acc}.fasta")
        if os.path.exists(p):
            return p
    return None


def db_id_for(acc, db_ids):
    """Database plasmid_id for an outbreak accession (versioned / RefSeq-prefixed forms)."""
    for cand in db_ids.get(acc, []):
        return cand
    return None


def build_release_order():
    """Training tree is built from the training table; reference plasmids are
    appended in accession order. Returns (training_tree_path, [(id, vector)])."""
    ref = pd.read_csv(REF_TSV, sep="\t", low_memory=False)
    z = np.load(REF_VECTORS, allow_pickle=True)
    vec = dict(zip(z["ids"], z["vectors"]))
    ref_ids = sorted(ref.loc[ref.source == "reference", "plasmid_id"])
    return [(pid, vec[pid]) for pid in ref_ids]


def training_vectors():
    out = {}
    for p in sorted(glob.glob(os.path.join(TRAINING_DIR, "*", "fastas", "*.fasta"))):
        for rec in SeqIO.parse(p, "fasta"):
            out.setdefault(rec.id, kmer_vector(str(rec.seq)))
    return out


def rebuild_without(exclude_id, ref_order, train_df, train_vecs):
    """Rebuild the release tree with one plasmid removed (same order otherwise)."""
    tree = FounderTree(PLIN_THRESHOLDS)
    codes = {}
    for pid in sorted(set(train_df.plasmid_id)):
        if pid != exclude_id:
            codes[pid] = tree.assign(train_vecs[pid], key=pid)[0]
    for pid, v in ref_order:
        if pid != exclude_id:
            codes[pid] = tree.assign(v, key=pid)[0]
    return tree, codes


def mates(codes, pid_code, level):
    """Database plasmids sharing the first `level` code positions."""
    pref = pid_code[:level]
    return frozenset(p for p, c in codes.items() if c[:level] == pref)


def main():
    meta = pd.read_csv(META_TSV, sep="\t")
    ref = pd.read_csv(REF_TSV, sep="\t", low_memory=False).drop_duplicates("plasmid_id")
    db_code = {pid: tuple(int(x) for x in r) for pid, r in zip(ref.plasmid_id, ref[BINS].values)}
    db_ids = {}
    for pid in ref.plasmid_id:
        base = pid.replace("RefSeq_", "").replace("NZ_", "").split(".")[0]
        db_ids.setdefault(base, []).append(pid)

    release = FounderTree.from_npz(REF_TREE)
    founders = set(release.founder_of.values())

    rows = []
    for _, m in meta.iterrows():
        acc = m["accession"]
        path = outbreak_fasta(acc)
        if path is None:
            print(f"  missing FASTA for {acc}; skipped")
            continue
        v = kmer_vector(str(next(SeqIO.parse(path, "fasta")).seq))
        code, _, new_level = release.assign(v, add=False)
        db_id = db_id_for(acc, db_ids)
        rows.append({
            "accession": acc, "study": m["study"], "resistance_gene": m["resistance_gene"],
            "country": m["country"], "source": m["source"], "predicted_inc": m["predicted_inc"],
            "pLIN": ".".join(map(str, code)), **dict(zip(BINS, code)),
            "in_database": db_id is not None, "database_id": db_id,
            "database_pLIN": ".".join(map(str, db_code[db_id])) if db_id else None,
            "is_new_code": new_level is not None,
            "is_founder": bool(db_id in founders) if db_id else None,
        })
    res = pd.DataFrame(rows)
    in_db = res[res.in_database]
    print(f"{len(res)} outbreak plasmids; {len(in_db)} already in the database")
    # the FASTA used here and the database copy should give the same code
    res["release_code_matches_database"] = np.where(
        res.in_database, res.pLIN == res.database_pLIN, None)

    # leave-one-out on database members; non-founders are unchanged by construction
    print("Leave-one-out rebuilds for founders ...")
    train_df = pd.read_csv(TRAIN_TSV, sep="\t")
    ref_order = build_release_order()
    tv = training_vectors()
    full_codes = {pid: c for pid, c in db_code.items()}
    for lvl in range(1, 7):
        res[f"loo_L{lvl}_concordant"] = None
    for i, r in res[res.in_database].iterrows():
        pid = r.database_id
        if not r.is_founder:
            for lvl in range(1, 7):
                res.at[i, f"loo_L{lvl}_concordant"] = True
            continue
        print(f"  rebuilding without {pid} ...")
        tree, codes = rebuild_without(pid, ref_order, train_df, tv)
        v = tv[pid] if pid in tv else dict(ref_order)[pid]
        loo_code, _, _ = tree.assign(v, add=False)
        others_full = {p: c for p, c in full_codes.items() if p != pid}
        for lvl in range(1, 7):
            res.at[i, f"loo_L{lvl}_concordant"] = (
                mates(codes, loo_code, lvl) == mates(others_full, full_codes[pid], lvl))

    # outbreak grouping
    pair_rows = []
    for (_, a), (_, b) in combinations(res.iterrows(), 2):
        pair_rows.append({"same_study": a.study == b.study,
                          **{f"share_L{l}": a.pLIN.split(".")[:l] == b.pLIN.split(".")[:l] for l in (3, 5, 6)}})
    pairs = pd.DataFrame(pair_rows)

    res.to_csv(OUT_TSV, sep="\t", index=False)
    loo = res[res.in_database]
    summary = {
        "n_plasmids": int(len(res)), "n_studies": int(res.study.nunique()),
        "n_in_database": int(len(loo)), "n_founders_rebuilt": int(loo.is_founder.sum()),
        "release_code_matches_database_copy": f"{int((loo.release_code_matches_database == True).sum())}/{len(loo)}",
        "loo_concordance": {f"L{l}": f"{int(loo[f'loo_L{l}_concordant'].astype(bool).sum())}/{len(loo)}"
                            for l in range(1, 7)},
        "n_new_codes_not_in_database": int(res.is_new_code.sum()),
        "grouping": {
            "same_study_pairs": int(pairs.same_study.sum()),
            "different_study_pairs": int((~pairs.same_study).sum()),
            **{f"pct_same_study_sharing_L{l}": float(100 * pairs.loc[pairs.same_study, f"share_L{l}"].mean())
               for l in (3, 5, 6)},
            **{f"pct_diff_study_sharing_L{l}": float(100 * pairs.loc[~pairs.same_study, f"share_L{l}"].mean())
               for l in (3, 5, 6)},
        },
    }
    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
