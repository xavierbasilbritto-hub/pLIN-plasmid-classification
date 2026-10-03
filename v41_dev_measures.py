#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
pLIN v4.1 development: which pairwise similarity separates related plasmids best?

DEVELOPMENT DATA ONLY: the 5,544 alignment-checked pairs of the v4 evaluation
(both halves) and the Swiss VIM-1 plasmids, all already seen. Results guide the
v4.1 design; they are not evidence of v4.1 performance (that comes from the
pre-registered confirmatory test on fresh plasmids).

Measures
  cosine            4-mer composition cosine distance (v3/v4 L5–L6)
  kmer_cont_short   k-mer (k=21, FracMinHash scaled 200) containment of the shorter
  kmer_min_cont     shared k-mers / larger sketch (symmetric; ~ bidirectional AF x identity^21)
  prot_cont         protein-family containment (v4 L1–L4)
  prot_jacc         protein-family Jaccard
  idf_cont          prevalence-weighted (log N/df) protein containment
  idf_jacc          prevalence-weighted protein Jaccard
Weighted ROC AUC (pair weights = inverse sampling probability) per truth.

Usage:
  python v41_dev_measures.py
Output: output/backbone_v41/dev/{pair_measures.tsv, auc.tsv}
"""

import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.metrics import roc_auc_score

import plin_v4 as v
from plin_kmers import containment, min_containment, sketch

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
V4 = os.path.join(BASE_DIR, "output", "backbone_v4")
DEV = os.path.join(BASE_DIR, "output", "backbone_v41", "dev")


def _sk(path):
    return sketch("".join(str(r.seq) for r in SeqIO.parse(path, "fasta")).upper())


def idf_weights():
    df = np.load(os.path.join(V4, "family_prevalence.npy"))
    n = 127517
    w = np.zeros(len(df) + 1)
    w[:len(df)] = np.log(n / np.maximum(df, 1))
    return w


def main():
    os.makedirs(DEV, exist_ok=True)
    pairs = pd.read_csv(os.path.join(V4, "truth_pairs.tsv"), sep="\t")
    ev = pd.read_csv(os.path.join(V4, "eval_plasmids.tsv"), sep="\t").set_index("plasmid_id")
    ids = sorted(set(pairs.plasmid_shorter) | set(pairs.plasmid_longer))
    with ProcessPoolExecutor(12) as ex:
        sk = dict(zip(ids, ex.map(_sk, [ev.fasta[i] for i in ids], chunksize=16)))
    fams, vecs = v.load_family_sets(), v.load_vectors()
    w = idf_weights()
    unit = lambda x: x / np.linalg.norm(x)

    rows = []
    for a, b in zip(pairs.plasmid_shorter, pairs.plasmid_longer):
        fa, fb = set(fams[a].tolist()), set(fams[b].tolist())
        inter, union = fa & fb, fa | fb
        wa, wb, wi, wu = (sum(w[x] for x in s) for s in (fa, fb, inter, union))
        rows.append({"cosine": 1 - float(unit(vecs[a]) @ unit(vecs[b])),
                     "kmer_cont_short": containment(sk[a], sk[b]),
                     "kmer_min_cont": min_containment(sk[a], sk[b]),
                     "prot_cont": len(inter) / min(len(fa), len(fb)) if fa and fb else 0.0,
                     "prot_jacc": len(inter) / len(union) if union else 0.0,
                     "idf_cont": wi / min(wa, wb) if wa and wb else 0.0,
                     "idf_jacc": wi / wu if wu else 0.0})
    m = pd.concat([pairs.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    m.to_csv(os.path.join(DEV, "pair_measures.tsv"), sep="\t", index=False)

    out = []
    for truth in ("same_lineage", "related_backbone"):
        for meas in ("cosine", "kmer_cont_short", "kmer_min_cont", "prot_cont", "prot_jacc", "idf_cont", "idf_jacc"):
            score = -m[meas] if meas == "cosine" else m[meas]
            out.append({"truth": truth, "measure": meas,
                        "AUC_weighted": roc_auc_score(m[truth], score, sample_weight=m.weight),
                        "AUC_unweighted": roc_auc_score(m[truth], score)})
    auc = pd.DataFrame(out)
    auc.to_csv(os.path.join(DEV, "auc.tsv"), sep="\t", index=False)
    print(auc.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
