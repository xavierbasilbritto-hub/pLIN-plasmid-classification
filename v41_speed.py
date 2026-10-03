#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
pLIN v4.1 typing speed (same machine and threads as benchmark_speed.py).

  known       100 release plasmids by whole-sequence match (instant path)
  related     the same 100 typed from sequence (proteins + k-mers; most proteins
              match the catalogue exactly)
  divergent   the 100 mutation-set plasmids with 1% random substitutions (seed 4343),
              so most proteins are unseen and go through the MMseqs2 search
The release is loaded once beforehand (load time reported separately); the
MMseqs2 search index is already built.

Usage:
  python v41_speed.py [--index-dir DIR]
Output: output/backbone_v41/speed_v41.json
"""

import argparse
import json
import os
import random
import time

import pandas as pd
from Bio import SeqIO

from plin_v41_typer import PlinV41Release
from v41_confirm_checks import mutate

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V41 = os.path.join(BASE_DIR, "output", "backbone_v41")
THREADS = 8


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index-dir", default=None)
    args = ap.parse_args()
    t = time.perf_counter()
    rel = PlinV41Release(os.path.join(V41, "release"), index_dir=args.index_dir)
    load_s = time.perf_counter() - t
    req = pd.read_csv(os.path.join(V41, "confirm", "requery_plasmids.tsv"), sep="\t").head(100)
    recs = [(i, "".join(str(r.seq) for r in SeqIO.parse(f, "fasta"))) for i, f in zip(req.plasmid_id, req.fasta)]
    mut = pd.read_csv(os.path.join(V41, "confirm", "mutation_plasmids.tsv"), sep="\t")
    r = random.Random(4343)
    div = [(f"{i}__1pct", mutate("".join(str(x.seq) for x in SeqIO.parse(f, "fasta")).upper(), 0.01, r))
           for i, f in zip(mut.plasmid_id, mut.fasta)]
    res = {"threads": THREADS, "release_load_seconds": round(load_s, 1)}
    for name, records, whole in (("known", recs, True), ("related", recs, False), ("divergent", div, False)):
        t = time.perf_counter()
        rel.type_sequences(records, threads=THREADS, whole_match=whole)
        s = time.perf_counter() - t
        res[name] = {"plasmids": len(records), "seconds": round(s, 2), "seconds_per_plasmid": round(s / len(records), 3)}
        print(name, res[name], flush=True)
    sp = pd.read_csv(os.path.join(BASE_DIR, "output", "backbone_v4", "speed", "speed_results.tsv"), sep="\t")
    mob = sp[(sp.tool == "MOB-suite") & (sp.source == "measured")]
    res["MOB-suite_seconds_per_plasmid_8_threads"] = round(float(mob.wall_s.sum() / mob.n_plasmids.sum()), 3)
    json.dump(res, open(os.path.join(V41, "speed_v41.json"), "w"), indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
