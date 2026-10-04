#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Counts for every step of the pLIN v4.1 database build-up (Supplementary Figure S1).

Sources -> merged IDs -> unique accessions -> proteins -> protein families ->
k-mer sketches -> founder tree -> release files. Every count is computed from the
build outputs; nothing is typed in by hand.

Usage:
  python v41_build_stats.py
Output: output/backbone_v41/build_stats.json
"""

import json
import os

import numpy as np
import pandas as pd

import plin_v4 as v

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V4 = os.path.join(BASE_DIR, "output", "backbone_v4")
V41 = os.path.join(BASE_DIR, "output", "backbone_v41")


def main():
    s = {}
    train = pd.read_csv(os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv"), sep="\t", usecols=["plasmid_id"])
    ref = pd.read_csv(os.path.join(BASE_DIR, "output", "pLIN_reference_assignments.tsv"), sep="\t", usecols=["plasmid_id"])
    s["training_rows"] = len(train)
    s["training_unique_ids"] = int(train.plasmid_id.nunique())
    s["reference_table_unique_ids"] = int(ref.plasmid_id.nunique())
    acc = ref.plasmid_id.drop_duplicates().str.replace("^RefSeq_", "", regex=True)
    s["unique_accessions"] = int(acc.nunique())
    s["accessions_listed_twice"] = int(s["reference_table_unique_ids"] - s["unique_accessions"])
    plsdb = set(pd.read_csv(os.path.join(V4, "onehealth", "plsdb_meta", "nuccore.csv"), usecols=["NUCCORE_ACC"]).NUCCORE_ACC)
    tr_acc = set(train.plasmid_id.str.replace("^RefSeq_", "", regex=True))
    ua = set(acc)
    s["source_PLSDB_2024_05_31"] = len(ua & plsdb)
    s["source_training_not_in_PLSDB"] = len((ua & tr_acc) - plsdb)
    s["source_other_NCBI"] = len(ua - plsdb - tr_acc)

    from v41_sketch_all import unique_ids
    ids = unique_ids()
    fs = v.load_family_sets()
    sizes = np.array([len(fs[i]) for i in ids])
    s["protein_families_per_plasmid_median"] = float(np.median(sizes))
    s["plasmids_without_proteins"] = int((sizes == 0).sum())
    s["proteins_predicted_all_ids"] = int(sum(1 for line in open(os.path.join(V4, "all_proteins.faa")) if line.startswith(">")))
    z = np.load(os.path.join(V4, "family_exact_lookup.npz") if os.path.exists(os.path.join(V4, "family_exact_lookup.npz"))
                else os.path.join(V4, "release", "family_exact_lookup.npz"))
    s["unique_protein_sequences"] = int(len(z["hash"]))
    s["protein_families_in_use"] = int(len(np.unique(z["family"])))
    fam = np.load(os.path.join(V4, "family_sets.npz"), allow_pickle=True)
    s["family_representatives"] = int(len(fam["family_rep"]))

    sk = np.load(os.path.join(V41, "sketches.npz"))
    lens = np.diff(sk["offsets"])
    s["kmer_hashes_total"] = int(lens.sum())
    s["kmer_hashes_per_plasmid_median"] = int(np.median(lens))

    dv = json.load(open(os.path.join(V41, "release", "DATABASE_VERSION.json")))
    s["release_unique_plasmids"] = dv["unique_plasmids"]
    s["clusters_per_level"] = dv["clusters_per_level"]
    s["release_files_bytes"] = dv["bytes"]
    hx = np.load(os.path.join(V41, "release", "plin_v41_index.npz"), allow_pickle=False)
    lv = hx["f_level"]
    s["founders_per_level"] = {f"L{k + 1}": int((lv == k).sum()) for k in range(4)}
    json.dump(s, open(os.path.join(V41, "build_stats.json"), "w"), indent=1)
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main()
