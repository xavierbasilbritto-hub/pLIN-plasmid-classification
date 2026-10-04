#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Build pLIN v4.1 codes for the release database (PREREGISTRATION_v4.1.md, section 2).

Every unique release plasmid (accession order) is coded with
plin_v41_nn.HybridIndex using the pre-registered levels:
  L1 prot 0.40, L2 prot 0.60, L3 kcont 0.50, L4 kcont 0.80  (earliest founder)
  L5 kmin 0.80, L6 kmin 0.95                               (nearest-neighbour copy)
Inputs: canonical protein-family sets (plin_v4.load_family_sets) and adaptive
k-mer sketches (output/backbone_v41/sketches.npz).

Output: output/backbone_v41/codes_v41.tsv, output/backbone_v41/hybrid_index.pkl
(working copy of the built index for the confirmatory checks).
"""

import os
import pickle
import time

import numpy as np
import pandas as pd

import plin_v4 as v
from plin_v41_nn import HybridIndex

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "output", "backbone_v41")
FOUNDER_LEVELS = [("prot", 0.40), ("prot", 0.60), ("kcont", 0.50), ("kcont", 0.80)]
LINEAGE = [0.80, 0.95]


def load_inputs():
    z = np.load(os.path.join(OUT, "sketches.npz"))
    ids, off, h, sc = z["ids"].tolist(), z["offsets"], z["hashes"], z["scaled"]
    S = [h[off[k]:off[k + 1]] for k in range(len(ids))]
    C = [int(x) for x in sc]
    fs = v.load_family_sets()
    F = [fs[i] for i in ids]
    return ids, F, S, C


def main():
    ids, F, S, C = load_inputs()
    t = time.perf_counter()
    hx = HybridIndex(FOUNDER_LEVELS, LINEAGE, F, S, C)
    hx.build(F, S, C, progress=20000)
    print(f"built {len(ids):,} codes in {time.perf_counter() - t:.0f} s", flush=True)
    codes = [".".join(map(str, c)) for c in hx.codes]
    pd.DataFrame({"plasmid_id": ids, "accession": [i.replace("RefSeq_", "", 1) for i in ids], "pLIN_v41": codes}) \
        .to_csv(os.path.join(OUT, "codes_v41.tsv"), sep="\t", index=False)
    with open(os.path.join(OUT, "hybrid_index.pkl"), "wb") as fh:
        pickle.dump({"ids": ids, "index": hx}, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print("clusters per level:", [len({c.rsplit('.', 5 - k)[0] if k < 5 else c for c in codes}) for k in range(6)])


if __name__ == "__main__":
    main()
