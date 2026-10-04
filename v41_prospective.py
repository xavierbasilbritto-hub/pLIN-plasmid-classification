#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Prospective surveillance simulation (secondary analysis, not pre-registered; development
data): the 74 published outbreak plasmids "arrive" in the order they were deposited in
NCBI, and each is typed when it arrives.

At arrival, the database holds only plasmids deposited on or before that day: every
release plasmid deposited later (and every plasmid identical to the arriving one) is
treated as absent, never copied from and never used as a founder. Isolates that
arrived earlier stay in the session, as in a hospital's running analysis.

For each isolate whose outbreak study already had an earlier isolate:
  early_warning   it shares the level with at least one earlier isolate of its study
  false_alarm     it shares the level with at least one earlier isolate of another study
MOB-suite clusters (fixed reference database, so independent of arrival order) are
scored on the same sequence of arrivals. Dates are NCBI deposition dates, not
sample-collection dates.

Usage:
  python v41_prospective.py [--index-dir DIR]
Output: output/backbone_v41/case_studies/{prospective_arrivals.tsv, prospective_metrics.tsv}
"""

import argparse
import glob
import os
import tempfile

import numpy as np
import pandas as pd
from Bio import SeqIO

from plin_kmers import adaptive_sketch, min_containment
from plin_v4 import seq_hash
from plin_v41_typer import PlinV41Release

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V41 = os.path.join(BASE_DIR, "output", "backbone_v41")
OUT = os.path.join(V41, "case_studies")
MOB_DIR = os.path.join(BASE_DIR, "output", "comparator_benchmark", "mob_typer")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index-dir", default=None)
    args = ap.parse_args()
    rel = PlinV41Release(os.path.join(V41, "release"), index_dir=args.index_dir)
    dates = pd.read_csv(os.path.join(V41, "accession_dates.tsv"), sep="\t", parse_dates=["create_date"])
    date_of = dict(zip(dates.accession, dates.create_date))
    db_dates = np.array([date_of.get(a, pd.NaT) for a in rel.accessions], dtype="datetime64[ns]")

    ob = pd.read_csv(os.path.join(BASE_DIR, "output", "outbreak_validation_founder_results.tsv"), sep="\t")
    files = {os.path.basename(f).rsplit(".", 1)[0]: f for f in
             glob.glob(os.path.join(BASE_DIR, "outbreak_validation", "*.fasta")) +
             glob.glob(os.path.join(BASE_DIR, "outbreak_validation", "expanded_sequences", "*.fasta"))}
    seqs = {a: "".join(str(r.seq) for r in SeqIO.parse(files[a], "fasta")).upper() for a in ob.accession}
    ob["date"] = [date_of.get(a, pd.NaT) for a in ob.accession]
    assert ob.date.notna().all(), "outbreak plasmid without a deposition date"
    ob = ob.sort_values(["date", "accession"]).reset_index(drop=True)

    row_of = {a: k for k, a in enumerate(rel.accessions)}
    h = pd.read_csv(os.path.join(rel.dir, "plasmid_hashes.tsv.gz"), sep="\t", dtype=str)
    rows_by_hash = h.groupby("hash").accession.apply(list).to_dict()
    fams = rel.families([(a, seqs[a]) for a in ob.accession], tempfile.mkdtemp(prefix="plin_prosp_"), threads=8)

    sess = rel.hx.session()
    s_sk, s_sc, s_codes, arrivals = [], [], [], []
    for _, r in ob.iterrows():
        a = r.accession
        later = set(np.flatnonzero(db_dates > np.datetime64(r.date)).tolist())   # not yet deposited
        same = {row_of[x] for x in rows_by_hash.get(str(np.uint64(seq_hash(seqs[a]))), []) if x in row_of}
        sk, sc = adaptive_sketch(seqs[a])
        extra = None
        if s_codes:
            extra = {"codes": np.array(s_codes), "kmin": np.array([min_containment(sk, x, sc, y) for x, y in zip(s_sk, s_sc)])}
        code, _ = sess.assign(fams[a], sk, sc, rel.n, extra=extra, add=True, exclude=later | same)
        s_sk.append(sk); s_sc.append(sc); s_codes.append(code)
        arrivals.append({"order": len(arrivals) + 1, "accession": a, "study": r.study, "date": r.date.date(),
                         "database_plasmids_available": int(rel.n - len(later)), "pLIN_v41_at_arrival": ".".join(map(str, code))})
    arr = pd.DataFrame(arrivals)

    clean = lambda x: None if pd.isna(x) or x in ("-", "") else x
    mob = {a: pd.read_csv(os.path.join(MOB_DIR, f"{a}.txt"), sep="\t", dtype=str).iloc[0] for a in arr.accession}
    arr["MOB_primary"] = [clean(mob[a].get("primary_cluster_id")) for a in arr.accession]
    arr["MOB_secondary"] = [clean(mob[a].get("secondary_cluster_id")) for a in arr.accession]
    arr.to_csv(os.path.join(OUT, "prospective_arrivals.tsv"), sep="\t", index=False)

    methods = {f"pLIN v4.1 L{k}": [".".join(c.split(".")[:k]) for c in arr.pLIN_v41_at_arrival] for k in (1, 3, 5, 6)}
    methods["MOB-suite primary"] = arr.MOB_primary.tolist()
    methods["MOB-suite secondary"] = arr.MOB_secondary.tolist()
    rows = []
    studies = arr.study.tolist()
    for name, lab in methods.items():
        warn, alarm, n_eval = [], [], 0
        for i in range(len(arr)):
            earlier_same = [j for j in range(i) if studies[j] == studies[i]]
            earlier_other = [j for j in range(i) if studies[j] != studies[i]]
            if not earlier_same:
                continue
            n_eval += 1
            warn.append(lab[i] is not None and any(lab[j] == lab[i] for j in earlier_same))
            alarm.append(lab[i] is not None and any(lab[j] == lab[i] for j in earlier_other))
        rows.append({"method": name, "isolates_with_earlier_outbreak_isolate": n_eval,
                     "early_warning_pct": round(100 * np.mean(warn), 1), "false_alarm_pct": round(100 * np.mean(alarm), 1)})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "prospective_metrics.tsv"), sep="\t", index=False)
    print(f"arrivals {arr.date.min()} to {arr.date.max()}; database available at first arrival: "
          f"{arr.database_plasmids_available.iloc[0]:,}, at last: {arr.database_plasmids_available.iloc[-1]:,}")
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
