#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Unit tests for the pLIN v4 founder tree (plin_v4.py).

  python test_plin_v4.py          # synthetic tests (seconds)
  python test_plin_v4.py --pilot  # also on the 2,747 protein-pilot plasmids (real families)
"""

import os
import sys
import unittest

import numpy as np
import pandas as pd

from plin_v4 import V4Tree, _ProteinNode, build_codes

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CHECK_PILOT = "--pilot" in sys.argv
L3 = 0.50


def synthetic(n=800, seed=0):
    """Family sets drawn from a few overlapping pools, 4-mer-like vectors around pool centres."""
    rng = np.random.default_rng(seed)
    pools = [rng.choice(3000, 120, replace=False) for _ in range(10)]
    centres = rng.dirichlet(np.ones(256) * 5, 10)
    fams, vecs = {}, {}
    for i in range(n):
        k = i % 10
        own = rng.choice(pools[k], int(rng.integers(5, 100)), replace=False)
        other = rng.choice(pools[(k + 1) % 10], int(rng.integers(0, 30)), replace=False)
        fams[f"P{i:05d}"] = np.unique(np.concatenate([own, other])).astype(np.int64)
        noise = rng.normal(0, rng.choice([0.0003, 0.003, 0.02]), 256) * centres[k]
        vecs[f"P{i:05d}"] = np.clip(centres[k] + noise, 1e-9, None)
    fams["P_EMPTY1"], vecs["P_EMPTY1"] = np.array([], dtype=np.int64), centres[0].copy()
    fams["P_EMPTY2"], vecs["P_EMPTY2"] = np.array([], dtype=np.int64), centres[0].copy()
    return fams, vecs


def containment(a, b):
    a, b = set(a.tolist()), set(b.tolist())
    return len(a & b) / min(len(a), len(b)) if a and b else 0.0


def cosine(u, v):
    return 1 - float(u @ v / np.linalg.norm(u) / np.linalg.norm(v))


def founders(ids, codes):
    """prefix -> founder id (the first plasmid in processing order to carry that prefix)."""
    f = {}
    for i in ids:
        for k in range(1, 7):
            f.setdefault(codes[i][:k], i)
    return f


class Checks:
    def check_all(self, fams, vecs, ids):
        codes, tree = build_codes(ids, fams, vecs, L3)

        # 1. re-querying any plasmid already in the tree returns its code
        for i in ids:
            self.assertEqual(tree.assign(fams[i], vecs[i], add=False), codes[i], i)

        # 2. every member is within its founder's threshold at every level
        fnd = founders(ids, codes)
        for i in ids:
            for li in range(6):
                f = fnd[codes[i][:li + 1]]
                if tree.kinds[li] == "protein":
                    if len(fams[i]) == 0:
                        self.assertEqual(codes[i][li], 0)
                    else:
                        self.assertGreaterEqual(containment(fams[i], fams[f]), tree.thresholds[li] - 1e-12)
                else:
                    self.assertLessEqual(cosine(vecs[i], vecs[f]), tree.thresholds[li] + 1e-12)

        # 3. adding plasmids never changes existing codes
        half = ids[: len(ids) // 2]
        codes_a, tree_a = build_codes(half, fams, vecs, L3)
        for i in ids[len(ids) // 2:]:
            tree_a.assign(fams[i], vecs[i])
        for i in half:
            self.assertEqual(tree_a.assign(fams[i], vecs[i], add=False), codes_a[i], i)

        # 4. same input, same order -> same codes
        self.assertEqual(build_codes(ids, fams, vecs, L3)[0], codes)
        return codes


class TestSynthetic(unittest.TestCase, Checks):
    def test_first_hit_matches_brute_force(self):
        rng = np.random.default_rng(1)
        node, sets = _ProteinNode(), []
        for c in range(1, 200):
            s = np.unique(rng.choice(400, int(rng.integers(1, 60)))).astype(np.int64)
            node.append(c, s)
            sets.append(s)
        for _ in range(500):
            q = np.unique(rng.choice(400, int(rng.integers(1, 60)))).astype(np.int64)
            t = float(rng.choice([0.2, 0.35, 0.5, 0.75]))
            expect = next((j for j, s in enumerate(sets) if containment(q, s) >= t - 1e-12), None)
            self.assertEqual(node.first_hit(q, t)[0], expect)

    def test_tree_guarantees(self):
        fams, vecs = synthetic()
        codes = self.check_all(fams, vecs, sorted(fams))
        self.assertEqual(codes["P_EMPTY1"][:4], (0, 0, 0, 0))
        self.assertEqual(codes["P_EMPTY1"], codes["P_EMPTY2"])     # identical vectors share the code


@unittest.skipUnless(CHECK_PILOT, "pass --pilot to run on the protein-pilot plasmids")
class TestPilot(unittest.TestCase, Checks):
    def test_pilot_guarantees(self):
        clu = pd.read_csv(os.path.join(BASE_DIR, "output", "protein_pilot", "clu_cluster.tsv"),
                          sep="\t", header=None, names=["rep", "member"])
        fam_id = {r: k for k, r in enumerate(sorted(clu.rep.unique()), 1)}
        clu["pid"] = clu.member.str.split("|").str[0]
        clu["f"] = clu.rep.map(fam_id)
        z = np.load(os.path.join(BASE_DIR, "output", "comparator_benchmark", "stability", "vectors.npz"),
                    allow_pickle=True)
        vecs = dict(zip(z["ids"].tolist(), z["X"].astype(np.float64)))
        fams = {p: np.unique(g.values).astype(np.int64) for p, g in clu.groupby("pid")["f"]}
        fams.update({p: np.array([], dtype=np.int64) for p in vecs if p not in fams})
        ids = sorted(vecs)
        codes = self.check_all(fams, vecs, ids)
        for k in range(1, 7):
            print(f"  L{k}: {len({c[:k] for c in codes.values()}):,} clusters for {len(ids):,} plasmids")


if __name__ == "__main__":
    unittest.main(argv=[a for a in sys.argv if a != "--pilot"], verbosity=2)
