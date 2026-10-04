#!/usr/bin/env python3
# Copyright (C) 2025-2026 Basil Britto Xavier. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Regenerate the inc_type column of the release with the current replicon classifier.

The column was inherited unchanged from the pLIN v4 release, so it predates the IncL/M group.
Until that group was added, pOXA-48-like plasmids were assigned to the nearest available type,
which was usually IncFII. This script reclassifies every release plasmid with the classifier in
data/inc_classifier.npz and writes a new dated release directory.

pLIN codes are not touched: every code, sketch, protein family and hash is copied unchanged, so
a code issued from the previous release stays valid. Only the replicon annotation changes.

Usage:
  python v41_reclassify_replicons.py [--release-dir DIR] [--out-dir DIR] [--version db-YYYY.MM.DD]
Output: a new release directory with plin_v41_codes.tsv.gz rebuilt, every other file copied,
        DATABASE_VERSION.json updated with new checksums, and reclassification_report.tsv
"""

import argparse
import hashlib
import json
import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.neighbors import KNeighborsClassifier

from build_inc_centroids import kmer_vector

BASE = os.path.dirname(os.path.abspath(__file__))
V41 = os.path.join(BASE, "output", "backbone_v41")
REF = os.path.join(BASE, "reference")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def vector_for(acc):
    """4-mer vector for one release accession, or None if its FASTA is not held locally."""
    p = os.path.join(REF, f"{acc}.fasta")
    if not os.path.exists(p):
        return acc, None
    try:
        seq = str(next(SeqIO.parse(p, "fasta")).seq)
    except Exception:
        return acc, None
    return acc, kmer_vector(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release-dir", default=os.path.join(V41, "release"))
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--version", default=None, help="new database version, e.g. db-2026.10.05")
    args = ap.parse_args()

    src = args.release_dir
    dv = json.load(open(os.path.join(src, "DATABASE_VERSION.json")))
    version = args.version or dv["database_version"]
    out = args.out_dir or os.path.join(V41, f"release_{version}")
    os.makedirs(out, exist_ok=True)

    clf = np.load(os.path.join(BASE, "data", "inc_classifier.npz"), allow_pickle=True)
    groups = list(clf["group_names"])
    knn = KNeighborsClassifier(n_neighbors=5, metric="cosine", weights="distance").fit(clf["X"], clf["y"])
    print(f"classifier: {len(groups)} groups, {clf['X'].shape[0]:,} training plasmids")

    codes = pd.read_csv(os.path.join(src, "plin_v41_codes.tsv.gz"), sep="\t", dtype=str)
    print(f"release: {len(codes):,} plasmids; reclassifying", flush=True)
    old = codes.inc_type.copy()

    accs = codes.accession.tolist()
    preds, missing = {}, 0
    with ProcessPoolExecutor(8) as ex:
        batch, keys = [], []
        for k, (acc, vec) in enumerate(ex.map(vector_for, accs, chunksize=64), 1):
            if vec is None:
                missing += 1
            else:
                batch.append(vec)
                keys.append(acc)
            if len(batch) >= 2000:
                for a, p in zip(keys, knn.predict(np.array(batch))):
                    preds[a] = groups[int(p)] if str(p).isdigit() else str(p)
                batch, keys = [], []
            if k % 20000 == 0:
                print(f"  {k:,}/{len(accs):,}", flush=True)
        if batch:
            for a, p in zip(keys, knn.predict(np.array(batch))):
                preds[a] = groups[int(p)] if str(p).isdigit() else str(p)

    codes["inc_type"] = codes.accession.map(preds).fillna(codes.inc_type)
    changed = (old != codes.inc_type) & old.notna()
    print(f"\nreclassified {len(preds):,} plasmids ({missing:,} without a local FASTA keep their previous label)")
    print(f"label changed for {int(changed.sum()):,} plasmids")
    rep = pd.DataFrame({"accession": codes.accession[changed], "previous": old[changed],
                        "current": codes.inc_type[changed]})
    rep.to_csv(os.path.join(out, "reclassification_report.tsv"), sep="\t", index=False)
    top = rep.groupby(["previous", "current"]).size().sort_values(ascending=False).head(8)
    print("\nlargest changes (previous -> current):")
    for (a, b), n in top.items():
        print(f"  {a} -> {b}: {n:,}")

    codes.to_csv(os.path.join(out, "plin_v41_codes.tsv.gz"), sep="\t", index=False, compression="gzip")
    for f in dv["sha256"]:
        if f != "plin_v41_codes.tsv.gz":
            shutil.copy(os.path.join(src, f), os.path.join(out, f))

    dv["database_version"] = version
    dv["supersedes"] = {"database_version": json.load(open(os.path.join(src, "DATABASE_VERSION.json")))["database_version"],
                        "change": "inc_type column regenerated with the IncL/M group added to the replicon "
                                  "classifier; pLIN codes, sketches, protein families and hashes are unchanged"}
    dv["replicon_classifier"] = {"groups": len(groups), "training_plasmids": int(clf["X"].shape[0]),
                                 "cv_accuracy": round(float(clf["cv_accuracy"].item()), 4)}
    dv["sha256"] = {f: sha256(os.path.join(out, f)) for f in dv["sha256"]}
    json.dump(dv, open(os.path.join(out, "DATABASE_VERSION.json"), "w"), indent=1)
    print(f"\nwrote {out}")
    print(f"new checksum for plin_v41_codes.tsv.gz: {dv['sha256']['plin_v41_codes.tsv.gz'][:16]}...")


if __name__ == "__main__":
    main()
