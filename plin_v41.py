# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
pLIN v4.1 founder tree: each level compares a plasmid with the founders of the
clusters inside its parent cluster, using one of

  prot   protein-family containment: shared families / families of the plasmid
         with fewer families, counted only when at least min(MIN_SHARED, that
         plasmid's family count) families are shared (one shared transposase no longer makes a tiny plasmid
         "contained" in anything)
  kcont  k-mer containment of the plasmid with fewer k-mers in the other
         (FracMinHash, plin_kmers; compared at the coarser of the two scales)
  kmin   shared k-mers / k-mers of the larger plasmid (symmetric: high only when
         both plasmids are mostly shared)

and joins either the earliest-created founder above the threshold ("earliest",
as in v3/v4) or the most similar one ("best"). Codes are nested; new clusters
get the next unused ID of their level, so existing codes never change.
"""

import numpy as np

from plin_kmers import MAX_SCALED, _cutoff

MIN_SHARED = 3
_SCALES = [2 ** i for i in range(int(np.log2(MAX_SCALED)) + 1)]          # 1, 2, ..., 256
_CUTS = np.array([_cutoff(s) for s in _SCALES], dtype=np.uint64)
_SCALE_IDX = {s: i for i, s in enumerate(_SCALES)}


def _pick(sim, t, rule):
    ok = np.flatnonzero(sim >= t - 1e-12)
    if not len(ok):
        return None
    return int(ok[0]) if rule == "earliest" else int(ok[np.argmax(sim[ok])])   # argmax: first on ties


class ProtNode:
    __slots__ = ("ids", "sizes", "post", "fp")

    def __init__(self):
        self.ids, self.sizes, self.post, self.fp = [], [], {}, []

    def similarity(self, fams):
        if not self.ids or not len(fams):
            return np.zeros(len(self.ids))
        hits = [self.post[f] for f in fams.tolist() if f in self.post]
        if not hits:
            return np.zeros(len(self.ids))
        inter = np.bincount(np.concatenate(hits), minlength=len(self.ids))
        small = np.minimum(np.asarray(self.sizes), len(fams))
        sim = inter / small
        sim[inter < np.minimum(MIN_SHARED, small)] = 0.0
        return sim

    def copy(self):
        n = ProtNode()
        n.ids, n.sizes, n.fp = list(self.ids), list(self.sizes), list(self.fp)
        n.post = {f: list(v) for f, v in self.post.items()}
        return n

    def append(self, cid, fams, fp=-1):
        j = len(self.ids)
        self.ids.append(cid)
        self.fp.append(fp)
        self.sizes.append(len(fams))
        for f in fams.tolist():
            self.post.setdefault(f, []).append(j)


class KmerNode:
    """Founder sketches with a lazily rebuilt sorted index (hash -> founder position)."""
    __slots__ = ("ids", "sk", "scale", "counts", "ix_hash", "ix_owner", "n_indexed", "fp")

    def __init__(self):
        self.ids, self.sk, self.scale, self.counts, self.fp = [], [], [], [], []
        self.ix_hash = np.zeros(0, np.uint64)
        self.ix_owner = np.zeros(0, np.int64)
        self.n_indexed = 0

    def copy(self):
        n = KmerNode()
        n.ids, n.sk, n.scale, n.counts, n.fp = list(self.ids), list(self.sk), list(self.scale), list(self.counts), list(self.fp)
        n.ix_hash, n.ix_owner, n.n_indexed = self.ix_hash, self.ix_owner, self.n_indexed   # arrays are replaced, never mutated
        return n

    def append(self, cid, sk, scale, fp=-1):
        self.ids.append(cid)
        self.fp.append(fp)
        self.sk.append(sk)
        self.scale.append(scale)
        self.counts.append(np.searchsorted(sk, _CUTS, side="right"))     # sketch size at every scale
        n = len(self.ids)
        if n - self.n_indexed > max(8, self.n_indexed // 4):
            h = np.concatenate(self.sk)
            o = np.repeat(np.arange(n), [s.size for s in self.sk])
            order = np.argsort(h, kind="stable")
            self.ix_hash, self.ix_owner, self.n_indexed = h[order], o[order], n

    def similarity(self, sk, scale, kind):
        n = len(self.ids)
        if not n or not sk.size:
            return np.zeros(n)
        inter = np.zeros(n)
        if self.n_indexed:
            lo = np.searchsorted(self.ix_hash, sk, side="left")
            hi = np.searchsorted(self.ix_hash, sk, side="right")
            m = hi > lo
            if m.any():
                idx = np.concatenate([np.arange(a, b) for a, b in zip(lo[m], hi[m])])
                inter += np.bincount(self.ix_owner[idx], minlength=n)[:n]
        for j in range(self.n_indexed, n):
            inter[j] = np.intersect1d(sk, self.sk[j], assume_unique=True).size
        # compare at the coarser scale: shared hashes are already below both cut-offs
        fs = np.asarray(self.scale)
        s_ix = np.array([_SCALE_IDX[max(scale, f)] for f in fs])
        q_cnt = np.searchsorted(sk, _CUTS, side="right")[s_ix]
        f_cnt = np.asarray(self.counts)[np.arange(n), s_ix]
        if kind == "kcont":
            den = np.minimum(q_cnt, f_cnt)
        else:                                                              # kmin
            den = np.maximum(q_cnt, f_cnt)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(den > 0, inter / den, 0.0)


class V41Tree:
    """levels: list of (kind, threshold); rule: 'earliest' or 'best'."""

    def __init__(self, levels, rule="best"):
        self.kinds = [k for k, _ in levels]
        self.thresholds = [t for _, t in levels]
        self.rule = rule
        self.children = {}
        self.next_id = [1] * len(levels)
        self._owned = None                 # session copies: prefixes already copied

    def session_copy(self):
        """Copy-on-write view: a session can add clusters without touching this tree."""
        c = V41Tree.__new__(V41Tree)
        c.kinds, c.thresholds, c.rule = list(self.kinds), list(self.thresholds), self.rule
        c.children, c.next_id, c._owned = dict(self.children), list(self.next_id), set()
        return c

    def _node_for_append(self, prefix, kind):
        node = self.children.get(prefix)
        if node is None:
            node = self.children[prefix] = ProtNode() if kind == "prot" else KmerNode()
            if self._owned is not None:
                self._owned.add(prefix)
        elif self._owned is not None and prefix not in self._owned:
            node = self.children[prefix] = node.copy()
            self._owned.add(prefix)
        return node

    def assign(self, fams, sk, scale, add=True, exclude=None):
        """exclude: database plasmid indices whose founder entries are ignored (leave-one-out)."""
        prefix, provisional = (), list(self.next_id)
        for li, (kind, t) in enumerate(zip(self.kinds, self.thresholds)):
            if kind == "prot" and len(fams) == 0:
                prefix += (0,)                     # no proteins: protein levels bucket 0
                continue
            node = self.children.get(prefix)
            chosen = None
            if node is not None:
                sim = node.similarity(fams) if kind == "prot" else node.similarity(sk, scale, kind)
                if exclude and node.fp:
                    sim = np.where(np.isin(np.asarray(node.fp), list(exclude)), -1.0, sim)
                j = _pick(sim, t, self.rule)
                chosen = node.ids[j] if j is not None else None
            if chosen is None:
                if add:
                    chosen = self.next_id[li]
                    self.next_id[li] += 1
                    node = self._node_for_append(prefix, kind)
                    if kind == "prot":
                        node.append(chosen, fams)
                    else:
                        node.append(chosen, sk, scale)
                else:
                    chosen = provisional[li]
                    provisional[li] += 1
            prefix += (chosen,)
        return prefix
