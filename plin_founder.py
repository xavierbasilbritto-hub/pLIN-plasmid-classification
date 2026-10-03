# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Founder-based (incremental, LIN-style) pLIN code assignment.

Why not single linkage: single linkage joins two plasmids at a level whenever
any chain of members links them, so clusters grow by chaining (at L5 one
cluster held 83.7% of the 8,077 training plasmids; at L3, 99.7%), and one new
"bridge" plasmid can merge two existing lineages. Why not average/complete
linkage: the tree is rebuilt on every run, so adding plasmids can split or
reshuffle existing clusters and earlier codes change.

Founder rule (leader clustering, as in CD-HIT): plasmids are processed in a
fixed order. At each level, within the cluster already chosen at the level
above, a plasmid joins the earliest-created cluster whose founder (first
member) is within the level's cosine-distance threshold; otherwise it founds
a new cluster with the next unused ID for that level. Consequences:
  * existing codes never change when plasmids are added (stability by design);
  * every member is within the threshold of its founder, so any two members
    of a cluster are within 2 x threshold of each other (no chaining);
  * the same input in the same order always gives the same codes;
  * any plasmid already in the tree re-queries to exactly its recorded code
    (this is why the earliest qualifying founder is used, not the nearest:
    a later founder could be nearer, which would make lookups disagree with
    the published table).

The state (founder vectors per cluster, next free ID per level) is saved
alongside the assignments so later releases and query-mode lookups extend the
same tree instead of rebuilding it.
"""

import numpy as np

LEVELS = ["A", "B", "C", "D", "E", "F"]

# pLIN hierarchical thresholds (cosine distance on 4-mer frequencies), the
# single source of truth for every pipeline (training, reference, app).
PLIN_THRESHOLDS = {
    "A": 0.150,   # L1 (Broad plasmid family)
    "B": 0.100,   # L2 (Subfamily)
    "C": 0.050,   # L3 (Cluster)
    "D": 0.020,   # L4 (Subcluster)
    "E": 0.010,   # L5 (Clone group)
    "F": 0.001,   # L6 (Lineage / Outbreak)
}

_BASE_CODE = np.full(256, -1, dtype=np.int64)
for _i, _b in enumerate(b"ACGT"):
    _BASE_CODE[_b] = _i
    _BASE_CODE[ord(chr(_b).lower())] = _i


def kmer_vector(sequence, k=4):
    """Tetranucleotide frequency vector: overlapping (sliding-window) counts of
    every ACGT-only 4-mer, divided by the number of such windows.

    This is the one vector definition shared by the classifier, the training
    and reference pipelines and the app (earlier versions of assign_pLIN.py and
    the app used str.count, which skips overlapping matches such as AAAA within
    AAAAA, so their vectors differed slightly from the classifier's).
    """
    codes = _BASE_CODE[np.frombuffer(sequence.encode("ascii", "replace"), dtype=np.uint8)]
    n = len(codes) - k + 1
    if n <= 0:
        return np.zeros(4 ** k)
    idx = np.zeros(n, dtype=np.int64)
    valid = np.ones(n, dtype=bool)
    for j in range(k):
        c = codes[j:j + n]
        valid &= c >= 0
        idx = idx * 4 + np.where(c >= 0, c, 0)
    counts = np.bincount(idx[valid], minlength=4 ** k).astype(np.float64)
    total = counts.sum()
    return counts / total if total > 0 else counts


def _unit(X):
    """Unit vector, rounded to float32 precision.

    The tree is stored as float32, so rounding here at build time makes every
    later reload compute bit-identical distances (otherwise a plasmid sitting
    exactly at a threshold could re-query differently after a save/load).
    """
    X = np.asarray(X, dtype=np.float64)
    n = np.linalg.norm(X, axis=-1, keepdims=True)
    n[n == 0] = 1.0
    return (X / n).astype(np.float32).astype(np.float64)


class _Node:
    """Founders of one cluster's children, in creation order (growable buffer)."""
    __slots__ = ("ids", "buf", "n")

    def __init__(self, ids=None, vecs=None):
        self.ids = list(ids) if ids is not None else []
        self.n = len(self.ids)
        cap = max(8, self.n)
        self.buf = np.zeros((cap, 256), dtype=np.float64)
        if self.n:
            self.buf[:self.n] = vecs

    def vecs(self):
        return self.buf[:self.n]

    def append(self, cid, u):
        if self.n == len(self.buf):
            grown = np.zeros((2 * len(self.buf), self.buf.shape[1]), dtype=np.float64)
            grown[:self.n] = self.buf[:self.n]
            self.buf = grown
        self.buf[self.n] = u
        self.n += 1
        self.ids.append(cid)

    def copy(self):
        return _Node(self.ids, self.vecs())


class FounderTree:
    """Hierarchical founder index: one node per cluster, keyed by its code prefix."""

    def __init__(self, thresholds):
        self.levels = list(thresholds.keys())
        self.thresholds = [float(thresholds[l]) for l in self.levels]
        self.children = {}                  # prefix tuple -> _Node
        self.next_id = [1] * len(self.levels)
        self.founder_of = {}                # full prefix tuple -> founder key
        self._owned = None                  # session copies: prefixes already copied

    def assign(self, vector, key=None, add=True):
        """Return (code tuple, distances to the chosen founder per level, new_level).

        new_level is the index of the first level at which a new cluster was
        founded (None if the plasmid joined existing clusters at every level).
        With add=False the tree is not modified and new clusters get provisional
        IDs (next_id + k) that are not reserved.
        """
        u = _unit(vector)
        prefix = ()
        dists, new_level = [], None
        provisional = list(self.next_id)
        for li, t in enumerate(self.thresholds):
            node = self.children.get(prefix)
            chosen, d_chosen = None, None
            if node is not None and node.n:
                # row-wise products: deterministic regardless of node size or
                # BLAS threading, so a build and a later lookup agree exactly
                d = 1.0 - (node.vecs() * u).sum(axis=1)
                hits = np.flatnonzero(d <= t)
                if len(hits):
                    # first founder in creation order, not the nearest: founders
                    # created later are always scanned later, so every plasmid
                    # already in the tree re-queries to exactly its recorded code
                    j = int(hits[0])
                    chosen, d_chosen = node.ids[j], float(d[j])
            if chosen is None:
                if new_level is None:
                    new_level = li
                if add:
                    chosen = self.next_id[li]
                    self.next_id[li] += 1
                    self._writable(prefix).append(chosen, u)
                    self.founder_of[prefix + (chosen,)] = key
                    d_chosen = 0.0
                else:
                    chosen = provisional[li]
                    provisional[li] += 1
            dists.append(d_chosen)
            prefix = prefix + (chosen,)
        return prefix, dists, new_level

    def _writable(self, prefix):
        node = self.children.get(prefix)
        if node is None:
            node = self.children[prefix] = _Node()
            if self._owned is not None:
                self._owned.add(prefix)
        elif self._owned is not None and prefix not in self._owned:
            node = self.children[prefix] = node.copy()      # copy-on-write
            self._owned.add(prefix)
        return node

    # ── persistence ──────────────────────────────────────────────────────────
    def to_npz(self, path, **extra):
        prefixes, child_ids, offsets, vec_rows = [], [], [0], []
        for prefix, node in self.children.items():
            prefixes.append(".".join(map(str, prefix)))
            child_ids.extend(node.ids)
            vec_rows.append(node.vecs())
            offsets.append(offsets[-1] + node.n)
        founders = sorted(self.founder_of.items())
        np.savez_compressed(
            path,
            levels=np.array(self.levels), thresholds=np.array(self.thresholds),
            next_id=np.array(self.next_id), prefixes=np.array(prefixes, dtype=object),
            child_ids=np.array(child_ids), offsets=np.array(offsets),
            vectors=np.vstack(vec_rows).astype(np.float32) if vec_rows else np.zeros((0, 256), np.float32),
            founder_prefix=np.array([".".join(map(str, p)) for p, _ in founders], dtype=object),
            founder_key=np.array([k for _, k in founders], dtype=object),
            **extra)

    @classmethod
    def from_npz(cls, path):
        z = np.load(path, allow_pickle=True)
        tree = cls(dict(zip(z["levels"].tolist(), z["thresholds"].tolist())))
        tree.next_id = z["next_id"].tolist()
        offsets, child_ids, vectors = z["offsets"], z["child_ids"], z["vectors"].astype(np.float64)
        for i, p in enumerate(z["prefixes"]):
            prefix = tuple(int(x) for x in p.split(".")) if p else ()
            s, e = offsets[i], offsets[i + 1]
            tree.children[prefix] = _Node([int(c) for c in child_ids[s:e]], vectors[s:e])
        for p, k in zip(z["founder_prefix"], z["founder_key"]):
            tree.founder_of[tuple(int(x) for x in p.split("."))] = k
        return tree


def canonical_order(keys):
    """Documented processing order for an initial build: sort by plasmid accession."""
    return sorted(range(len(keys)), key=lambda i: keys[i])


def build_founder_codes(vectors, keys, thresholds, order=None):
    """Assign founder codes to all plasmids; returns (codes list of tuples, tree)."""
    tree = FounderTree(thresholds)
    order = canonical_order(keys) if order is None else order
    codes = [None] * len(keys)
    for i in order:
        codes[i], _, _ = tree.assign(vectors[i], key=keys[i])
    return codes, tree


def session_copy(tree):
    """Copy a tree so a batch of queries can extend it without touching the release tree.

    Nodes are shared until a query appends to one, which then gets its own
    copy (copy-on-write), so concurrent sessions never see each other's queries.
    """
    c = FounderTree.__new__(FounderTree)
    c.levels, c.thresholds = list(tree.levels), list(tree.thresholds)
    c.children = dict(tree.children)
    c.next_id = list(tree.next_id)
    c.founder_of = dict(tree.founder_of)
    c._owned = set()
    return c
