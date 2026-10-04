#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Query-side inputs for v4.1 development (protein families via the release's
nearest-catalogued-protein search, and adaptive k-mer sketches), cached once:

  vim        19 Swiss VIM-1 plasmid contigs
  outbreak   74 published outbreak plasmids
  mutants    50 release plasmids (seed 3) with random point mutations at
             0.01%, 0.1% and 1% (seed 7), plus the unmutated originals

DEVELOPMENT DATA ONLY (all already used during v4 work).
Output: output/backbone_v41/dev/query_inputs.npz
"""

import glob
import os
import random

import numpy as np
import pandas as pd
from Bio import SeqIO

import plin_v4 as v
from plin_kmers import adaptive_sketch
from plin_v4_typer import PlinV4Release

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DEV = os.path.join(BASE_DIR, "output", "backbone_v41", "dev")


def mutate(s, rate, r):
    s = list(s)
    for k in r.sample(range(len(s)), max(1, int(len(s) * rate))):
        s[k] = r.choice([b for b in "ACGT" if b != s[k]])
    return "".join(s)


def main():
    os.makedirs(DEV, exist_ok=True)
    recs = []
    for f in sorted(glob.glob(os.path.join(BASE_DIR, "sample_data", "swiss_vim1_outbreak", "*_plasmids.fasta"))):
        recs += [("vim", r.id, str(r.seq).upper()) for r in SeqIO.parse(f, "fasta")]
    ob = pd.read_csv(os.path.join(BASE_DIR, "output", "outbreak_validation_founder_results.tsv"), sep="\t")
    files = {os.path.basename(f).rsplit(".", 1)[0]: f for f in
             glob.glob(os.path.join(BASE_DIR, "outbreak_validation", "*.fasta")) +
             glob.glob(os.path.join(BASE_DIR, "outbreak_validation", "expanded_sequences", "*.fasta"))}
    recs += [("outbreak", a, "".join(str(r.seq) for r in SeqIO.parse(files[a], "fasta")).upper()) for a in ob.accession]
    p = v.release_plasmids()
    ids = random.Random(3).sample([i for i in sorted(p) if not i.startswith("RefSeq_")], 50)
    orig = {i: str(next(SeqIO.parse(p[i], "fasta")).seq).upper() for i in ids}
    recs += [("mut_orig", i, orig[i]) for i in ids]
    for rate in (0.0001, 0.001, 0.01):
        r = random.Random(7)
        recs += [(f"mut_{rate}", f"{i}__{rate}", mutate(orig[i], rate, r)) for i in ids]

    rel = PlinV4Release(os.path.join(BASE_DIR, "output", "backbone_v4", "release"))
    work = os.path.join(DEV, "typer_work")
    os.makedirs(work, exist_ok=True)
    paths = {}
    for k, (_, rid, seq) in enumerate(recs):
        path = os.path.join(work, f"q{k}.fa")
        open(path, "w").write(f">{rid}\n{seq}\n")
        paths[rid] = path
    fams = rel._families(paths, work, threads=12)
    sks = [adaptive_sketch(seq) for _, _, seq in recs]
    lens = [len(fams[rid]) for _, rid, _ in recs]
    sk_lens = [s.size for s, _ in sks]
    np.savez(os.path.join(DEV, "query_inputs.npz"),
             group=np.array([g for g, _, _ in recs]), rid=np.array([r for _, r, _ in recs]),
             fam_off=np.concatenate([[0], np.cumsum(lens)]), fam=np.concatenate([fams[r] for _, r, _ in recs]),
             sk_off=np.concatenate([[0], np.cumsum(sk_lens)]), sk=np.concatenate([s for s, _ in sks]),
             scale=np.array([s for _, s in sks]), length=np.array([len(seq) for _, _, seq in recs]))
    print(pd.Series([g for g, _, _ in recs]).value_counts().to_dict())


if __name__ == "__main__":
    main()
