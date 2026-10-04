#!/usr/bin/env python3
# Copyright (C) 2025-2026 Basil Britto Xavier. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Which similarity measure separates plasmid relationships best? (secondary analysis, not pre-registered)

An earlier version of pLIN grouped plasmids by cosine distance on tetranucleotide (4-mer) frequency
vectors. v4.1 uses protein-family containment and FracMinHash k-mer sketches instead. This scores
the measures themselves, rather than whole schemes, on the confirmatory pairs and the registered
alignment truth, so the choice of measure can be judged independently of thresholds and of the
assignment rule.

Measures compared, for each aligned pair:
  cosine      cosine similarity of 256-dimensional 4-mer frequency vectors (the earlier method)
  prot_cont   shared protein families / families of the plasmid with fewer (v4.1 L1-L2)
  prot_jacc   shared protein families / families of either plasmid
  kcont       k-mer containment of the smaller sketch in the larger (v4.1 L3-L4)
  kmin        shared hashes / larger sketch, symmetric (v4.1 L5-L6)

Reported as the area under the ROC curve against each registered truth, which is threshold-free.

Usage:
  python v41_measure_comparison.py
Output: output/backbone_v41/confirm/measure_comparison.tsv and .json
"""

import json
import os

import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.metrics import roc_auc_score

from build_inc_centroids import kmer_vector
from plin_kmers import adaptive_sketch, containment, min_containment

BASE = os.path.dirname(os.path.abspath(__file__))
V41 = os.path.join(BASE, "output", "backbone_v41")
C = os.path.join(V41, "confirm")
REF = os.path.join(BASE, "reference")

MEASURES = [("cosine", "4-mer cosine similarity (earlier method)"),
            ("prot_cont", "protein-family containment"),
            ("prot_jacc", "protein-family Jaccard"),
            ("kcont", "k-mer containment"),
            ("kmin", "symmetric k-mer similarity")]


def sequence_measures(pairs):
    """cosine, kcont and kmin for every pair whose FASTA files are held locally."""
    ids = sorted(set(pairs.plasmid_shorter) | set(pairs.plasmid_longer))
    data = {}
    for n, i in enumerate(ids, 1):
        p = os.path.join(REF, f"{i}.fasta")
        if not os.path.exists(p):
            continue
        try:
            seq = str(next(SeqIO.parse(p, "fasta")).seq)
            h, sc = adaptive_sketch(seq)
            data[i] = (kmer_vector(seq), h, sc)
        except Exception:
            continue
        if n % 400 == 0:
            print(f"  {n}/{len(ids)}", flush=True)
    rows = []
    for a, b in zip(pairs.plasmid_shorter, pairs.plasmid_longer):
        if a not in data or b not in data:
            rows.append((np.nan, np.nan, np.nan))
            continue
        va, ha, sa = data[a]
        vb, hb, sb = data[b]
        cos = float(np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb)))
        rows.append((cos, containment(ha, hb, sa, sb), min_containment(ha, hb, sa, sb)))
    return pd.DataFrame(rows, columns=["cosine", "kcont", "kmin"])


def protein_measures(pairs):
    """Protein-family containment and Jaccard from the family sets of the database build."""
    fs = np.load(os.path.join(BASE, "output", "backbone_v4", "family_sets.npz"), allow_pickle=True)
    idx = {a: k for k, a in enumerate(list(fs["ids"]))}
    off, fam = fs["offsets"], fs["families"]

    def families(a):
        k = idx.get(a)
        return set(fam[off[k]:off[k + 1]].tolist()) if k is not None else None

    cont, jacc = [], []
    for a, b in zip(pairs.plasmid_shorter, pairs.plasmid_longer):
        A, B = families(a), families(b)
        if not A or not B:
            cont.append(np.nan)
            jacc.append(np.nan)
            continue
        shared = len(A & B)
        cont.append(shared / min(len(A), len(B)))
        jacc.append(shared / len(A | B))
    return pd.DataFrame({"prot_cont": cont, "prot_jacc": jacc})


def main():
    pairs = pd.read_csv(os.path.join(C, "truth_pairs.tsv"), sep="\t")
    print(f"{len(pairs):,} confirmatory pairs; computing measures", flush=True)
    d = pd.concat([pairs.reset_index(drop=True), sequence_measures(pairs), protein_measures(pairs)], axis=1)
    d = d.dropna(subset=[m for m, _ in MEASURES])
    d.to_csv(os.path.join(C, "measure_comparison.tsv"), sep="\t", index=False)

    out = {"pairs": int(len(d)), "same_lineage_pairs": int(d.same_lineage.sum()),
           "related_backbone_pairs": int(d.related_backbone.sum()), "measures": {}}
    print(f"\n{len(d):,} pairs with every measure\n")
    print(f'{"measure":42s} {"AUC lineage":>12s} {"AUC backbone":>13s}')
    print("-" * 70)
    for col, name in MEASURES:
        al = float(roc_auc_score(d.same_lineage, d[col]))
        ab = float(roc_auc_score(d.related_backbone, d[col]))
        out["measures"][col] = {"name": name, "AUC_same_lineage": round(al, 3),
                                "AUC_related_backbone": round(ab, 3)}
        print(f"{name:42s} {al:12.3f} {ab:13.3f}")
    json.dump(out, open(os.path.join(C, "measure_comparison.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
