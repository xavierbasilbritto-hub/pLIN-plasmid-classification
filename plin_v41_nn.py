# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
pLIN v4.1, nearest-neighbour (LIN-code style) assignment.

As in LIN codes for bacterial genomes (Hennart et al.), a plasmid copies the
code of its most similar earlier plasmid level by level: at level k, among the
plasmids that share its code at levels 1..k-1, the most similar one (ties: the
earliest) is found with that level's measure; if the similarity reaches the
level threshold its level-k ID is copied, otherwise a new ID is created and
every deeper level also gets a new ID. Existing codes never change, a plasmid
already in the database finds itself (similarity 1) and gets its own code back,
and a slightly changed plasmid copies the code of its original.

Measures (plin_v41): prot (protein-family containment, at least
min(MIN_SHARED, smaller set) families shared), kcont (k-mer containment of
the plasmid with fewer k-mers), kmin (shared k-mers / larger sketch).

Indexes are global inverted indexes over all indexed plasmids (protein family
-> plasmids; k-mer hash -> plasmids), so similarities to every candidate come
from one pass over the query's postings.
"""

import numpy as np

from plin_kmers import MAX_SCALED, _cutoff

MIN_SHARED = 3
_SCALES = [2 ** i for i in range(int(np.log2(MAX_SCALED)) + 1)]
_CUTS = np.array([_cutoff(s) for s in _SCALES], dtype=np.uint64)
_SIDX = {s: i for i, s in enumerate(_SCALES)}


class CSR:
    """Static inverted index key -> owner rows, from per-owner sorted key arrays."""

    def __init__(self, keys_per_owner):
        lens = np.array([len(k) for k in keys_per_owner])
        keys = np.concatenate(keys_per_owner) if len(keys_per_owner) else np.zeros(0)
        owners = np.repeat(np.arange(len(keys_per_owner), dtype=np.int32), lens)
        order = np.argsort(keys, kind="stable")
        self.keys, self.owners = keys[order], owners[order]

    def owners_of(self, q):
        lo = np.searchsorted(self.keys, q, side="left")
        hi = np.searchsorted(self.keys, q, side="right")
        m = hi > lo
        if not m.any():
            return np.zeros(0, np.int32)
        lo, lens = lo[m], (hi - lo)[m]
        start = np.repeat(lo, lens)                              # expand every [lo, hi) range
        off = np.arange(lens.sum()) - np.repeat(np.cumsum(lens) - lens, lens)
        return self.owners[start + off]


class NNIndex:
    """Database of plasmids (families, sketches, codes) with per-level nearest-neighbour coding."""

    def __init__(self, levels, fams, sketches, scales):
        self.kinds = [k for k, _ in levels]
        self.thr = np.array([t for _, t in levels])
        self.n = len(fams)
        self.fsize = np.array([len(f) for f in fams])
        self.prot = CSR([f.astype(np.int64) for f in fams])
        self.kmer = CSR([s for s in sketches])
        self.scale = np.asarray(scales)
        self.scale_ix = np.array([_SIDX[int(x)] for x in scales], dtype=np.int64)
        self.cnt = np.array([np.searchsorted(s, _CUTS, side="right") for s in sketches])   # size at each scale
        self.codes = np.zeros((self.n, len(levels)), dtype=np.int64)
        self.next_id = np.ones(len(levels), dtype=np.int64)

    def similarities(self, fams, sk, scale, limit):
        """Sparse similarities to database plasmids [0, limit): {kind: (sorted row indices, values)}.

        Only plasmids sharing at least one protein family / k-mer can pass a threshold > 0, so only
        those are scored.
        """
        out = {}
        if "prot" in self.kinds:
            rows, inter = np.unique(self.prot.owners_of(fams.astype(np.int64)), return_counts=True)
            keep = rows < limit
            rows, inter = rows[keep], inter[keep].astype(float)
            small = np.minimum(self.fsize[rows], len(fams))
            ok = (small > 0) & (inter >= np.minimum(MIN_SHARED, small))
            out["prot"] = (rows[ok], inter[ok] / small[ok])
        if "kcont" in self.kinds or "kmin" in self.kinds:
            rows, inter = np.unique(self.kmer.owners_of(sk), return_counts=True)
            keep = rows < limit
            rows, inter = rows[keep], inter[keep].astype(float)
            s_ix = np.maximum(self.scale_ix[rows], _SIDX[scale])
            qc = np.searchsorted(sk, _CUTS, side="right")[s_ix]
            dc = self.cnt[rows, s_ix]
            lo, hi = np.minimum(qc, dc), np.maximum(qc, dc)
            out["kcont"] = (rows, np.where(lo > 0, inter / np.maximum(lo, 1), 0.0))
            out["kmin"] = (rows, np.where(hi > 0, inter / np.maximum(hi, 1), 0.0))
        return out

    def code_for(self, sims, limit, noprot, extra=None, next_id=None):
        """Per-level nearest-neighbour code against database plasmids [0, limit).

        extra: optional session plasmids {"codes": (m, L) array, "sims": {kind: (m,) array}} that the
        query may also copy from (plasmids typed earlier in the same session).
        Returns (code list, index of the first newly created level or None).
        """
        next_id = self.next_id if next_id is None else next_id
        m = len(extra["codes"]) if extra else 0
        out, new_from = [], None
        for k, kind in enumerate(self.kinds):
            if new_from is None and kind == "prot" and noprot:
                out.append(0)                                   # no proteins: protein levels bucket 0
                continue
            if new_from is None:
                rows, vals = sims[kind]
                if k:
                    same = (self.codes[rows, :k] == np.asarray(out)).all(axis=1)
                    rows, vals = rows[same], vals[same]
                j = int(np.argmax(vals)) if len(vals) else -1   # rows ascending: first max = earliest
                best = float(vals[j]) if len(vals) else -1.0
                ebest, ecode = -1.0, None
                if m:
                    es = extra["sims"][kind].copy()
                    if k:
                        es[~(extra["codes"][:, :k] == np.asarray(out)).all(axis=1)] = -1.0
                    ej = int(np.argmax(es))
                    ebest, ecode = float(es[ej]), int(extra["codes"][ej, k])
                if max(best, ebest) >= self.thr[k] - 1e-12:
                    out.append(int(self.codes[rows[j], k]) if best >= ebest else ecode)
                    continue
                new_from = k
            out.append(int(next_id[k]))
            next_id[k] += 1
        return out, new_from

    def build(self, fams, sketches, scales, progress=0):
        """Code the indexed plasmids in index order (each against the earlier ones)."""
        for i in range(self.n):
            sims = self.similarities(fams[i], sketches[i], scales[i], i)
            code, _ = self.code_for(sims, i, noprot=len(fams[i]) == 0)
            self.codes[i] = code
            if progress and (i + 1) % progress == 0:
                print(f"  {i + 1:,}/{self.n:,}", flush=True)


class HybridIndex:
    """Founder rule for the backbone levels, nearest-neighbour copying for the lineage levels.

    levels: [(kind, t), ...] for the founder levels (L1..Lf), then lineage thresholds lin = [t5, t6, ...]
    on the symmetric k-mer similarity (kmin). A plasmid whose nearest earlier plasmid (kmin, ties:
    earliest) reaches lin[0] copies that plasmid's whole backbone code and its first lineage ID, then
    copies the next lineage ID from the nearest plasmid inside that lineage cluster while kmin reaches
    the next threshold. A plasmid with no such relative is placed in the backbone levels by the
    founder rule (plin_v41.V41Tree, earliest founder) and starts new lineage IDs.
    """

    def __init__(self, founder_levels, lin, fams, sketches, scales):
        from plin_v41 import V41Tree
        self.nf, self.lin = len(founder_levels), list(lin)
        self.tree = V41Tree(founder_levels, rule="earliest")
        self.ix = NNIndex([("kmin", 0.0)], [np.zeros(0, np.int64)] * len(fams), sketches, scales)
        self.ix.kinds = ["kmin"]
        self.codes = np.zeros((len(fams), self.nf + len(lin)), dtype=np.int64)
        self.next_lin = np.ones(len(lin), dtype=np.int64)

    def _lineage(self, rows, vals, codes_of):
        """Copy lineage IDs level by level from the nearest plasmid inside the current cluster."""
        out = []
        for k, t in enumerate(self.lin):
            if k:
                keep = codes_of[:, self.nf + k - 1] == out[-1]
                rows, vals, codes_of = rows[keep], vals[keep], codes_of[keep]
            if len(vals) and vals.max() >= t - 1e-12:
                out.append(int(codes_of[int(np.argmax(vals)), self.nf + k]))
            else:
                return out, k
        return out, None

    def assign(self, fams, sk, scale, limit, extra=None, add=True):
        """Code against indexed plasmids [0, limit) (+ session extras). Returns (code, new_from_level)."""
        rows, vals = self.ix.similarities(np.zeros(0, np.int64), sk, scale, limit)["kmin"]
        codes_of = self.codes[rows]
        if extra is not None and len(extra["codes"]):
            rows = np.concatenate([rows, -1 - np.arange(len(extra["codes"]))])
            vals = np.concatenate([vals, extra["kmin"]])
            codes_of = np.vstack([codes_of, extra["codes"]])
        next_lin = self.next_lin if add else self.next_lin.copy()
        if len(vals) and vals.max() >= self.lin[0] - 1e-12:
            j = int(np.argmax(vals))
            backbone = [int(x) for x in codes_of[j, :self.nf]]
            same_bb = (codes_of[:, :self.nf] == backbone).all(axis=1)
            lin, stop = self._lineage(rows[same_bb], vals[same_bb], codes_of[same_bb])
            new_from = None if stop is None else self.nf + stop
        else:
            before = list(self.tree.next_id)
            backbone = list(self.tree.assign(fams, sk, scale, add=add))
            lin = []
            new_from = next((k for k, c in enumerate(backbone) if c >= before[k]), self.nf)
        for k in range(len(lin), len(self.lin)):
            lin.append(int(next_lin[k]))
            next_lin[k] += 1
        return backbone + lin, new_from

    def session(self):
        """A typing session: shares the release index, but new clusters live only in the session."""
        import copy
        s = copy.copy(self)
        s.tree = self.tree.session_copy()
        s.next_lin = self.next_lin.copy()
        return s

    def build(self, fams, sketches, scales, progress=0):
        for i in range(len(fams)):
            code, _ = self.assign(fams[i], sketches[i], scales[i], i)
            self.codes[i] = code
            if progress and (i + 1) % progress == 0:
                print(f"  {i + 1:,}/{len(fams):,}", flush=True)


# ── compact release format ───────────────────────────────────────────────────

def save_hybrid(path, hx, fams, sketches, scales):
    """Release file: database codes and sketches (lineage search) + founder replay for L1..Lf.

    The founder of a backbone cluster is the first plasmid (index order) carrying its code prefix;
    loading re-appends founders in creation order (ascending cluster ID per node), which repeats the
    build's arithmetic exactly.
    """
    nf = hx.nf
    first = {}
    for i, c in enumerate(hx.codes):
        c = tuple(int(x) for x in c[:nf])
        for k in range(nf):
            if hx.tree.kinds[k] == "prot" and c[k] == 0:
                continue                                       # no-protein bucket: not a cluster
            first.setdefault(c[:k + 1], i)
    rows = sorted(first.items(), key=lambda kv: (len(kv[0]), kv[0][-1]))
    prot = [i for p, i in rows if hx.tree.kinds[len(p) - 1] == "prot"]
    kmer = [i for p, i in rows if hx.tree.kinds[len(p) - 1] != "prot"]
    sk_len = np.array([len(s) for s in sketches])
    np.savez(path,
             founder_levels=np.array([f"{k}:{t}" for k, t in zip(hx.tree.kinds, hx.tree.thresholds)]),
             lineage=np.array(hx.lin), codes=hx.codes, next_id=np.array(hx.tree.next_id), next_lin=hx.next_lin,
             f_level=np.array([len(p) - 1 for p, _ in rows], dtype=np.int8),
             f_parent=np.array([".".join(map(str, p[:-1])) for p, _ in rows]),
             f_cid=np.array([p[-1] for p, _ in rows], dtype=np.int64),
             f_plasmid=np.array([i for _, i in rows], dtype=np.int64),
             prot_off=np.concatenate([[0], np.cumsum([len(fams[i]) for i in prot])]).astype(np.int64),
             prot_val=np.concatenate([fams[i] for i in prot]).astype(np.int64),
             sk_off=np.concatenate([[0], np.cumsum(sk_len)]).astype(np.int64),
             sk_val=np.concatenate(sketches).astype(np.uint64), sk_scale=np.asarray(scales, dtype=np.int64))


def load_hybrid(path):
    from plin_v41 import V41Tree
    z = np.load(path, allow_pickle=False)
    levels = [(s.split(":")[0], float(s.split(":")[1])) for s in z["founder_levels"].tolist()]
    off, val, scales = z["sk_off"], z["sk_val"], z["sk_scale"]
    sketches = [val[off[k]:off[k + 1]] for k in range(len(scales))]
    hx = HybridIndex.__new__(HybridIndex)
    hx.nf, hx.lin = len(levels), [float(x) for x in z["lineage"]]
    hx.tree = V41Tree(levels, rule="earliest")
    hx.ix = NNIndex([("kmin", 0.0)], [np.zeros(0, np.int64)] * len(scales), sketches, scales)
    hx.ix.kinds = ["kmin"]
    hx.codes = z["codes"]
    hx.next_lin = z["next_lin"].copy()
    po, pv = z["prot_off"], z["prot_val"]
    pi = 0
    for lv, par, cid, pl in zip(z["f_level"], z["f_parent"], z["f_cid"], z["f_plasmid"]):
        prefix = tuple(int(x) for x in str(par).split(".")) if par else ()
        kind = levels[lv][0]
        node = hx.tree._node_for_append(prefix, kind)
        if kind == "prot":
            node.append(int(cid), pv[po[pi]:po[pi + 1]])
            pi += 1
        else:
            node.append(int(cid), sketches[pl], int(scales[pl]))
    hx.tree.next_id = z["next_id"].tolist()
    return hx
