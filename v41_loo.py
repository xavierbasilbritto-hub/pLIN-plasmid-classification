#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Leave-one-out outbreak validation of pLIN v4.1 (secondary analysis, not pre-registered;
mirrors the leave-one-out test of the earlier manuscript).

The 74 published outbreak plasmids (27 studies) are each typed from sequence,
one at a time, against the release with that plasmid and every database
plasmid with an identical sequence treated as absent (never copied from and
never used as a founder): the situation of a hospital receiving a new isolate.
The outbreak set was used during the design of v4.1, so this is a development
analysis.

Per level L1..L6:
  recovery          share of the database plasmids that get back their published code
  same_study_linked share of plasmids sharing the level with at least one other
                    plasmid of their own outbreak study (other plasmids keep their
                    published code, or their typed code if not in the database)
  false_link        share of plasmids sharing the level with any plasmid of another study
                    (several studies describe the same internationally spread plasmid, so
                    some of these links are real; the pairwise rates below are the
                    measure comparable with the outbreak benchmark)
  same_study_pairs  share of (held-out, other) same-study pairs sharing the level
  diff_study_pairs  share of (held-out, other) different-study pairs sharing the level

Usage:
  python v41_loo.py [--index-dir DIR]
Output: output/backbone_v41/case_studies/{loo_codes.tsv, loo_metrics.tsv}
"""

import argparse
import glob
import os
import tempfile

import numpy as np
import pandas as pd
from Bio import SeqIO

from plin_kmers import adaptive_sketch
from plin_v4 import seq_hash
from plin_v41_typer import PlinV41Release

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "output", "backbone_v41", "case_studies")


def shared(a, b):
    n = 0
    for x, y in zip(str(a).split("."), str(b).split(".")):
        if x != y:
            break
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index-dir", default=None)
    args = ap.parse_args()
    rel = PlinV41Release(os.path.join(BASE_DIR, "output", "backbone_v41", "release"), index_dir=args.index_dir)
    row_of = {a: k for k, a in enumerate(rel.accessions)}
    h = pd.read_csv(os.path.join(rel.dir, "plasmid_hashes.tsv.gz"), sep="\t", dtype=str)
    rows_by_hash = h.groupby("hash").accession.apply(list).to_dict()

    ob = pd.read_csv(os.path.join(BASE_DIR, "output", "outbreak_validation_founder_results.tsv"), sep="\t")
    files = {os.path.basename(f).rsplit(".", 1)[0]: f for f in
             glob.glob(os.path.join(BASE_DIR, "outbreak_validation", "*.fasta")) +
             glob.glob(os.path.join(BASE_DIR, "outbreak_validation", "expanded_sequences", "*.fasta"))}
    seqs = {a: "".join(str(r.seq) for r in SeqIO.parse(files[a], "fasta")).upper() for a in ob.accession}
    study = dict(zip(ob.accession, ob.study))

    typed = pd.read_csv(os.path.join(OUT, "outbreak_codes.tsv"), sep="\t").set_index("plasmid_id")
    reference = {}                                   # code each plasmid has without leave-one-out
    excl, in_db = {}, {}
    for a in ob.accession:
        same = rows_by_hash.get(str(np.uint64(seq_hash(seqs[a]))), [])
        excl[a] = {row_of[x] for x in same if x in row_of}
        in_db[a] = bool(excl[a])
        reference[a] = typed.at[a, "pLIN_v41"]

    work = tempfile.mkdtemp(prefix="plin_loo_")
    fams = rel.families([(a, seqs[a]) for a in ob.accession], work, threads=8)
    loo = {}
    for a in ob.accession:
        sk, sc = adaptive_sketch(seqs[a])
        code, _ = rel.hx.session().assign(fams[a], sk, sc, rel.n, add=False, exclude=excl[a])
        loo[a] = ".".join(map(str, code))
    pd.DataFrame({"accession": ob.accession, "study": ob.study, "in_database": [in_db[a] for a in ob.accession],
                  "excluded_database_rows": [len(excl[a]) for a in ob.accession],
                  "published_or_typed_code": [reference[a] for a in ob.accession],
                  "loo_code": [loo[a] for a in ob.accession]}) \
        .to_csv(os.path.join(OUT, "loo_codes.tsv"), sep="\t", index=False)

    rows = []
    dbs = [a for a in ob.accession if in_db[a]]
    for k in range(1, 7):
        rec = np.mean([shared(loo[a], reference[a]) >= k for a in dbs])
        has_partner = [a for a in ob.accession if any(study[o] == study[a] for o in ob.accession if o != a)]
        linked = np.mean([any(shared(loo[a], reference[o]) >= k for o in ob.accession if o != a and study[o] == study[a])
                          for a in has_partner])
        false = np.mean([any(shared(loo[a], reference[o]) >= k for o in ob.accession if study[o] != study[a])
                         for a in ob.accession])
        same_p = [shared(loo[a], reference[o]) >= k for a in ob.accession for o in ob.accession
                  if o != a and study[o] == study[a]]
        diff_p = [shared(loo[a], reference[o]) >= k for a in ob.accession for o in ob.accession if study[o] != study[a]]
        rows.append({"level": f"L{k}", "plasmids_in_database": len(dbs), "recovery_pct": round(100 * rec, 1),
                     "plasmids_with_outbreak_partner": len(has_partner),
                     "same_study_linked_pct": round(100 * linked, 1), "false_link_pct": round(100 * false, 1),
                     "same_study_pairs_pct": round(100 * np.mean(same_p), 1),
                     "diff_study_pairs_pct": round(100 * np.mean(diff_p), 2)})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "loo_metrics.tsv"), sep="\t", index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
