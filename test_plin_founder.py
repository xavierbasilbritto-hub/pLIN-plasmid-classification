#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Unit tests for the guarantees pLIN codes rely on (plin_founder.py).

  python test_plin_founder.py            # synthetic tests (seconds)
  python test_plin_founder.py --release  # also checks the shipped release files
"""

import os
import sys
import tempfile
import unittest

import numpy as np

from plin_founder import (PLIN_THRESHOLDS, FounderTree, canonical_order, kmer_vector,
                          session_copy)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CHECK_RELEASE = "--release" in sys.argv


def synthetic_vectors(n=600, n_centres=12, seed=0):
    """4-mer-like vectors around a few centres, with spreads spanning the thresholds."""
    rng = np.random.default_rng(seed)
    centres = rng.dirichlet(np.ones(256) * 5, n_centres)
    X = []
    for i in range(n):
        c = centres[i % n_centres]
        noise = rng.normal(0, rng.choice([0.0005, 0.003, 0.01, 0.03]), 256) * c
        X.append(np.clip(c + noise, 1e-9, None))
    keys = [f"P{i:05d}" for i in range(n)]
    return np.array(X), keys


def build(X, keys, order=None):
    tree = FounderTree(PLIN_THRESHOLDS)
    codes = {}
    for i in (canonical_order(keys) if order is None else order):
        codes[keys[i]] = tree.assign(X[i], key=keys[i])[0]
    return tree, codes


def cos(u, v):
    return 1 - float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)))


class FounderGuarantees(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.X, cls.keys = synthetic_vectors()

    def test_deterministic(self):
        _, a = build(self.X, self.keys)
        _, b = build(self.X, self.keys)
        self.assertEqual(a, b)

    def test_requery_returns_recorded_code(self):
        tree, codes = build(self.X, self.keys)
        for i, k in enumerate(self.keys):
            self.assertEqual(tree.assign(self.X[i], add=False)[0], codes[k], k)

    def test_requery_after_save_and_load(self):
        tree, codes = build(self.X, self.keys)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "tree.npz")
            tree.to_npz(path)
            loaded = FounderTree.from_npz(path)
        self.assertEqual(loaded.next_id, tree.next_id)
        for i, k in enumerate(self.keys):
            self.assertEqual(loaded.assign(self.X[i], add=False)[0], codes[k], k)

    def test_existing_codes_never_change_when_plasmids_are_added(self):
        half = len(self.keys) // 2
        tree, before = build(self.X[:half], self.keys[:half])
        for i in range(half, len(self.keys)):
            tree.assign(self.X[i], key=self.keys[i])
        for i, k in enumerate(self.keys[:half]):
            self.assertEqual(tree.assign(self.X[i], add=False)[0], before[k], k)

    def test_members_within_threshold_of_founder(self):
        tree, codes = build(self.X, self.keys)
        idx = {k: i for i, k in enumerate(self.keys)}
        for prefix, founder in tree.founder_of.items():
            level = len(prefix) - 1
            t = list(PLIN_THRESHOLDS.values())[level]
            f = self.X[idx[founder]]
            for k, c in codes.items():
                if c[:level + 1] == prefix:
                    self.assertLessEqual(cos(self.X[idx[k]], f), t + 1e-6)

    def test_codes_are_nested(self):
        _, codes = build(self.X, self.keys)
        parent = {}
        for c in codes.values():
            for lvl in range(1, 6):
                key = (lvl, c[lvl])
                self.assertEqual(parent.setdefault(key, c[:lvl]), c[:lvl])

    def test_session_copy_leaves_release_tree_untouched(self):
        tree, codes = build(self.X[:300], self.keys[:300])
        snapshot = (list(tree.next_id), {p: n.n for p, n in tree.children.items()})
        s = session_copy(tree)
        for i in range(300, 600):
            s.assign(self.X[i], key=self.keys[i])
        self.assertEqual(snapshot, (list(tree.next_id), {p: n.n for p, n in tree.children.items()}))
        for i, k in enumerate(self.keys[:300]):
            self.assertEqual(tree.assign(self.X[i], add=False)[0], codes[k])

    def test_kmer_vector_counts_overlapping_windows(self):
        v = kmer_vector("AAAAAC")          # windows AAAA, AAAA, AAAC
        self.assertAlmostEqual(v[0], 2 / 3)
        self.assertAlmostEqual(v[1], 1 / 3)
        self.assertAlmostEqual(kmer_vector("acgtNacgt").sum(), 1.0)   # N-containing windows skipped


@unittest.skipUnless(CHECK_RELEASE, "pass --release to check the shipped release files")
class ReleaseFiles(unittest.TestCase):
    def test_reference_tree_reproduces_published_codes(self):
        import pandas as pd
        tree = FounderTree.from_npz(os.path.join(BASE_DIR, "data", "plin_founder_tree_reference.npz"))
        ref = pd.read_csv(os.path.join(BASE_DIR, "output", "pLIN_reference_assignments.tsv"),
                          sep="\t", low_memory=False)
        z = np.load(os.path.join(BASE_DIR, "output", "reference_kmer_vectors.npz"), allow_pickle=True)
        vec = dict(zip(z["ids"], z["vectors"]))
        sample = ref[ref.source == "reference"].sample(2000, random_state=0)
        for pid, code in zip(sample.plasmid_id, sample.pLIN):
            got = ".".join(map(str, tree.assign(vec[pid], add=False)[0]))
            self.assertEqual(got, code, pid)

    def test_training_codes_preserved_in_release(self):
        import pandas as pd
        tr = pd.read_csv(os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv"), sep="\t")
        ref = pd.read_csv(os.path.join(BASE_DIR, "output", "pLIN_reference_assignments.tsv"),
                          sep="\t", low_memory=False)
        m = tr.merge(ref[ref.source == "training"][["plasmid_id", "inc_type", "pLIN"]],
                     on=["plasmid_id", "inc_type"], suffixes=("", "_release"))
        self.assertEqual(len(m), len(tr))
        self.assertTrue((m.pLIN == m.pLIN_release).all())


if __name__ == "__main__":
    unittest.main(argv=[a for a in sys.argv if a != "--release"], verbosity=2)
