#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""FracMinHash sketches (k=21, adaptive power-of-two scale; plin_kmers) for every unique release plasmid.

One row per accession (RefSeq_ duplicates of db-2026.10.02 dropped, keeping the
plain accession). Output: output/backbone_v41/sketches.npz (ids, offsets, hashes).
"""
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from Bio import SeqIO

import plin_v4 as v
from plin_kmers import adaptive_sketch

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "backbone_v41")


def unique_ids():
    ids = sorted(v.release_plasmids())
    seen, keep = set(), []
    for i in sorted(ids, key=lambda x: (x.startswith("RefSeq_"), x)):
        a = i.replace("RefSeq_", "", 1)
        if a not in seen:
            seen.add(a)
            keep.append(i)
    return sorted(keep)


def _one(args):
    pid, path = args
    return adaptive_sketch("".join(str(r.seq) for r in SeqIO.parse(path, "fasta")).upper())


def main():
    paths = v.release_plasmids()
    ids = unique_ids()
    with ProcessPoolExecutor(14) as ex:
        res = list(ex.map(_one, [(i, paths[i]) for i in ids], chunksize=64))
    sks, scales = [r[0] for r in res], np.array([r[1] for r in res])
    lens = np.array([s.size for s in sks])
    np.savez(os.path.join(OUT, "sketches.npz"), ids=np.array(ids), offsets=np.concatenate([[0], np.cumsum(lens)]),
             hashes=np.concatenate(sks).astype(np.uint64), scaled=scales)
    print(f"{len(ids):,} plasmids, {lens.sum():,} hashes (median {int(np.median(lens))} per plasmid, "
          f"{(lens < 25).sum():,} plasmids with < 25)")


if __name__ == "__main__":
    main()
