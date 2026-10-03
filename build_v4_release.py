#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Assemble the frozen pLIN v4 database release (output/backbone_v4/release/).

Contents
  plin_v4_codes.tsv.gz        one row per unique accession: v4 code and levels,
                              v3 code (v3 -> v4 mapping), replicon group,
                              length, number of protein families, aliases
                              (5,788 accessions were present twice in
                              db-2026.10.02, as X and RefSeq_X, with identical
                              sequence and code; they are merged here)
  plin_v4_tree.npz            compact founder tree (plin_v4.load_tree)
  family_representatives.faa.gz   one sequence per protein family (MMseqs2
                              search target; `plin_v4.py index` builds the
                              search index from it)
  family_rep_to_family.tsv.gz representative name -> canonical family ID
  family_exact_lookup.npz     protein-sequence hash -> family (deviation 3)
  plasmid_hashes.tsv.gz       whole-plasmid sequence hash -> accession
  DATABASE_VERSION.json       counts, thresholds, software versions, SHA-256

Every code in the table is checked against the compact tree before writing.

Usage:
  python build_v4_release.py
"""

import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import plin_v4 as v

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V4 = os.path.join(BASE_DIR, "output", "backbone_v4")
REL = os.path.join(V4, "release")
RELEASE = "db-2026.10.03"
SCHEME = "pLIN v4.0"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    os.makedirs(REL, exist_ok=True)

    # verify the compact tree against the code table
    fams, vecs = v.load_family_sets(), v.load_vectors()
    codes = pd.read_csv(os.path.join(V4, "codes", f"plin_v4_codes_L3_{v.RELEASE_L3:.2f}.tsv"), sep="\t")
    tree = v.load_tree(v.TREE_FILE)
    bad = [i for i, c in zip(codes.plasmid_id, codes.pLIN_v4)
           if ".".join(map(str, tree.assign(fams[i], vecs[i], add=False))) != c]
    assert not bad, f"{len(bad)} codes not reproduced by the compact tree, e.g. {bad[:3]}"

    # one row per accession
    v3 = pd.read_csv(os.path.join(BASE_DIR, "output", "pLIN_reference_assignments.tsv"), sep="\t",
                     usecols=["plasmid_id", "pLIN", "inc_type", "length_bp"]).drop_duplicates("plasmid_id")
    d = codes.merge(v3, on="plasmid_id", how="left")
    d["accession"] = d.plasmid_id.str.replace("^RefSeq_", "", regex=True)
    conflicts = d.groupby("accession").agg(n4=("pLIN_v4", "nunique"), n3=("pLIN", "nunique"))
    assert (conflicts.n4 == 1).all(), "duplicate accessions with different v4 codes"
    alias = d.groupby("accession").plasmid_id.apply(lambda s: ";".join(sorted(s)))
    d = d.sort_values(["accession", "plasmid_id"]).drop_duplicates("accession")
    d["ids_in_db_2026_10_02"] = d.accession.map(alias)
    lv = d.pLIN_v4.str.split(".", expand=True)
    for k in range(6):
        d[f"L{k + 1}"] = lv[k].astype(int)
    out = d[["accession", "pLIN_v4"] + [f"L{k + 1}" for k in range(6)] +
            ["n_families", "pLIN", "inc_type", "length_bp", "ids_in_db_2026_10_02"]] \
        .rename(columns={"pLIN": "pLIN_v3"})
    out.to_csv(os.path.join(REL, "plin_v4_codes.tsv.gz"), sep="\t", index=False, compression="gzip")

    shutil.copy(v.TREE_FILE, os.path.join(REL, "plin_v4_tree.npz"))
    with open(os.path.join(V4, "clu_rep_seq.fasta"), "rb") as src, \
            gzip.open(os.path.join(REL, "family_representatives.faa.gz"), "wb", compresslevel=6) as dst:
        shutil.copyfileobj(src, dst)
    z = np.load(os.path.join(V4, "family_sets.npz"), allow_pickle=True)
    pd.DataFrame({"representative": z["family_rep"], "family": z["family_rep_family"]}) \
        .to_csv(os.path.join(REL, "family_rep_to_family.tsv.gz"), sep="\t", index=False, compression="gzip")
    shutil.copy(v.EXACT_PROTEINS, os.path.join(REL, "family_exact_lookup.npz"))
    pd.read_csv(v.EXACT_PLASMIDS, sep="\t", dtype=str) \
        .to_csv(os.path.join(REL, "plasmid_hashes.tsv.gz"), sep="\t", index=False, compression="gzip")

    import pyrodigal
    mmseqs = subprocess.run([v.MMSEQS, "version"], capture_output=True, text=True).stdout.strip()
    files = sorted(f for f in os.listdir(REL) if f != "DATABASE_VERSION.json")
    info = {
        "scheme": SCHEME, "database_version": RELEASE,
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "unique_plasmids": int(len(out)), "accessions_merged_from_duplicates": int(len(codes) - len(out)),
        "plasmids_without_proteins": int((out.n_families == 0).sum()),
        "clusters_per_level": {f"L{k + 1}": int(out[f"L{k + 1}"].nunique()) for k in range(6)},
        "unique_codes": int(out.pLIN_v4.nunique()),
        "thresholds": {"L1-L4 protein-family containment (>=)": [0.20, 0.35, v.RELEASE_L3, 0.75],
                       "L5-L6 4-mer cosine distance (<=)": [0.010, 0.001]},
        "assignment": "founder rule (earliest founder within threshold, nested by level); codes never renumbered",
        "software": {"pyrodigal": pyrodigal.__version__, "mmseqs2": mmseqs, "python": platform.python_version(),
                     "numpy": np.__version__},
        "protein_families": {"clustering": "MMseqs2 easy-cluster --min-seq-id 0.5 -c 0.8 --cov-mode 0",
                             "identical_sequences": "one family per identical sequence (lowest ID)",
                             "n_representatives": int(len(z["family_rep"]))},
        "preregistration": "https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification/blob/"
                           "preregistration-v4/output/backbone_v4/PREREGISTRATION.md",
        "sha256": {f: sha256(os.path.join(REL, f)) for f in files},
        "bytes": {f: os.path.getsize(os.path.join(REL, f)) for f in files},
    }
    json.dump(info, open(os.path.join(REL, "DATABASE_VERSION.json"), "w"), indent=2)
    print(json.dumps({k: info[k] for k in ("unique_plasmids", "accessions_merged_from_duplicates",
                                            "plasmids_without_proteins", "clusters_per_level", "unique_codes",
                                            "bytes")}, indent=2))


if __name__ == "__main__":
    main()
