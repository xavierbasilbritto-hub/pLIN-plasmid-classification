#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Assemble and verify the pLIN v4.1 database release (output/backbone_v41/release/).

Contents
  plin_v41_codes.tsv.gz      one row per unique accession: v4.1 code and levels L1–L6,
                             pLIN v4 and v3 codes (mapping), replicon group, length,
                             protein-family count, IDs in db-2026.10.02
  plin_v41_index.npz         database codes, k-mer sketches (lineage search) and the
                             backbone founder tree (plin_v41_nn.save_hybrid)
  family_members.faa.gz      every unique catalogued protein (MMseqs2 search target;
                             header = sequence hash)
  family_exact_lookup.npz    protein-sequence hash -> family
  plasmid_hashes.tsv.gz      whole-plasmid sequence hash -> accession
  DATABASE_VERSION.json      counts, levels, software versions, SHA-256, links

R4 (pre-registered release criterion): the index file alone must reproduce every
database code, and no accession may appear twice; the build stops otherwise.

Usage:
  python build_v41_release.py
"""

import hashlib
import multiprocessing
import json
import os
import pickle
import platform
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import plin_v4 as v
from plin_v41_nn import load_hybrid, save_hybrid

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V41 = os.path.join(BASE_DIR, "output", "backbone_v41")
V4REL = os.path.join(BASE_DIR, "output", "backbone_v4", "release")
REL = os.path.join(V41, "release")
RELEASE, SCHEME = "db-2026.10.03", "pLIN v4.1"
PREREG = ("https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification/blob/"
          "preregistration-v4/output/backbone_v41/PREREGISTRATION_v4.1.md")

_G = {}


def _check(rng):
    hx, F, S, Cs, codes = _G["hx"], _G["F"], _G["S"], _G["C"], _G["codes"]
    n = len(F)
    return [i for i in range(*rng) if hx.assign(F[i], S[i], Cs[i], n, add=False)[0] != codes[i]]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    os.makedirs(REL, exist_ok=True)
    with open(os.path.join(V41, "hybrid_index.pkl"), "rb") as fh:
        d = pickle.load(fh)
    ids, hx = d["ids"], d["index"]
    acc = [i.replace("RefSeq_", "", 1) for i in ids]
    assert len(set(acc)) == len(acc), "R4: duplicate accessions in the release"
    z = np.load(os.path.join(V41, "sketches.npz"))
    assert z["ids"].tolist() == ids
    off, h, sc = z["offsets"], z["hashes"], z["scaled"]
    S = [h[off[k]:off[k + 1]] for k in range(len(ids))]
    C = [int(x) for x in sc]
    fs = v.load_family_sets()
    F = [fs[i] for i in ids]

    idx_path = os.path.join(REL, "plin_v41_index.npz")
    save_hybrid(idx_path, hx, F, S, C)

    # R4: the file alone reproduces every code
    _G.update(hx=load_hybrid(idx_path), F=F, S=S, C=C, codes=[list(map(int, c)) for c in hx.codes])
    chunks = [(a, min(a + 2000, len(ids))) for a in range(0, len(ids), 2000)]
    with ProcessPoolExecutor(14, mp_context=multiprocessing.get_context("fork")) as ex:   # workers share _G
        bad = [i for part in ex.map(_check, chunks) for i in part]
    assert not bad, f"R4: {len(bad)} codes not reproduced from the release file, e.g. {[ids[i] for i in bad[:3]]}"
    print(f"R4: release file reproduces all {len(ids):,} codes; no duplicate accessions", flush=True)

    codes = pd.DataFrame({"plasmid_id": ids, "accession": acc,
                          "pLIN_v41": [".".join(map(str, c)) for c in hx.codes]})
    for k in range(6):
        codes[f"L{k + 1}"] = hx.codes[:, k]
    v4 = pd.read_csv(os.path.join(V4REL, "plin_v4_codes.tsv.gz"), sep="\t", dtype=str).set_index("accession")
    codes["n_families"] = [len(f) for f in F]
    codes["pLIN_v4"] = codes.accession.map(v4.pLIN_v4)
    codes["pLIN_v3"] = codes.accession.map(v4.pLIN_v3)
    codes["inc_type"] = codes.accession.map(v4.inc_type)
    codes["length_bp"] = codes.accession.map(v4.length_bp)
    codes["ids_in_db_2026_10_02"] = codes.accession.map(v4.ids_in_db_2026_10_02)
    codes.drop(columns="plasmid_id").to_csv(os.path.join(REL, "plin_v41_codes.tsv.gz"), sep="\t", index=False,
                                            compression="gzip")
    for f in ("family_members.faa.gz", "family_exact_lookup.npz"):
        shutil.copy(os.path.join(V4REL, f), os.path.join(REL, f))
    ph = pd.read_csv(os.path.join(V4REL, "plasmid_hashes.tsv.gz"), sep="\t", dtype=str)
    ph["accession"] = ph.plasmid_id.str.replace("^RefSeq_", "", regex=True)
    ph.drop_duplicates("accession")[["accession", "hash"]] \
        .to_csv(os.path.join(REL, "plasmid_hashes.tsv.gz"), sep="\t", index=False, compression="gzip")

    import pyrodigal
    mmseqs = subprocess.run([v.MMSEQS, "version"], capture_output=True, text=True).stdout.strip()
    files = sorted(f for f in os.listdir(REL) if f != "DATABASE_VERSION.json" and not os.path.isdir(os.path.join(REL, f)))
    info = {
        "scheme": SCHEME, "database_version": RELEASE,
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "unique_plasmids": len(ids), "plasmids_without_proteins": int(sum(len(f) == 0 for f in F)),
        "clusters_per_level": {f"L{k + 1}": int(len({tuple(c[:k + 1]) for c in hx.codes.tolist()})) for k in range(6)},
        "levels": {"L1": "backbone family: protein-family containment >= 0.40 (founder)",
                   "L2": "backbone group: protein-family containment >= 0.60 (founder)",
                   "L3": "shared backbone: k-mer containment >= 0.50 (founder)",
                   "L4": "backbone variant: k-mer containment >= 0.80 (founder)",
                   "L5": "lineage: symmetric k-mer similarity >= 0.80 (nearest-neighbour copy)",
                   "L6": "near-identical / outbreak clone: symmetric k-mer similarity >= 0.95 (nearest-neighbour copy)"},
        "kmers": "canonical 21-mers, splitmix64, FracMinHash with power-of-two scale (>= 400 kept, max 256)",
        "software": {"pyrodigal": pyrodigal.__version__, "mmseqs2": mmseqs, "python": platform.python_version(),
                     "numpy": np.__version__},
        "release_checks": "R1–R4 passed; see PREREGISTRATION_v4.1 and output/backbone_v41/confirm/",
        "preregistration": PREREG,
        "sha256": {f: sha256(os.path.join(REL, f)) for f in files},
        "bytes": {f: os.path.getsize(os.path.join(REL, f)) for f in files},
    }
    json.dump(info, open(os.path.join(REL, "DATABASE_VERSION.json"), "w"), indent=2)
    print(json.dumps({k: info[k] for k in ("unique_plasmids", "plasmids_without_proteins", "clusters_per_level", "bytes")},
                     indent=2))


if __name__ == "__main__":
    main()
