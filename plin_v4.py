#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
pLIN scheme v4: protein-family coarse levels (L1–L4) + 4-mer fine levels (L5–L6).

Pre-registered in output/backbone_v4/PREREGISTRATION.md (GitHub branch
preregistration-v4, pushed 2026-10-02 14:37 UTC). Summary:

  L1–L4  a plasmid joins the earliest founder (in creation order) whose
         protein-family containment with it is >= the level's threshold;
         containment = shared families / families of the plasmid with fewer
         families. Thresholds L1 0.20, L2 0.35, L3 calibrated, L4 0.75.
  L5–L6  as in v3: earliest founder within cosine distance 0.010 / 0.001 on
         4-mer frequency vectors.
  Every level searches only inside the cluster chosen at the level above, so
  codes are nested, and the founder rule makes existing codes permanent.

Plasmids with no predicted protein cannot be compared by containment; they
get L1–L4 = 0 and are coded at L5–L6 by 4-mer composition within that bucket
(a case the pre-registration did not cover; recorded under Deviations).

Steps (each resumable):
  predict     pyrodigal (metagenomic mode) per plasmid + 4-mer vectors. The
              Prodigal binary used in the pilot crashes on ~7% of plasmids
              (deviation 2 in the pre-registration).
  catalogue   MMseqs2 easy-cluster, >= 50% identity, >= 80% coverage
              (--cov-mode 0, as in the pilot); frozen family IDs
  canonicalize  one family per identical protein sequence (deviation 3)
  build       founder codes for the whole release database for one L3 value
  index       searchable representatives + release tree (query mode)
  query       codes for new plasmids (exact plasmid match, exact protein
              match, then MMseqs2 search)

Usage:
  python plin_v4.py predict --threads 12
  python plin_v4.py catalogue --threads 14
  python plin_v4.py build --l3 0.50
"""

import argparse
import glob
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from Bio import SeqIO

from plin_founder import _Node, _unit, kmer_vector
from validate_alignment_backbone import fasta_paths

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "output", "backbone_v4")
FAA = os.path.join(OUT, "faa_pyrodigal")
MMSEQS = os.path.expanduser("~/miniconda3/envs/pLIN_tools/bin/mmseqs")

PROTEIN_THRESHOLDS = {"A": 0.20, "B": 0.35, "C": None, "D": 0.75}   # containment, >=
KMER_THRESHOLDS = {"E": 0.010, "F": 0.001}                          # cosine distance, <=
L3_GRID = [0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70]


# ── inputs ────────────────────────────────────────────────────────────────────

def release_plasmids():
    """plasmid_id -> FASTA path for every plasmid in the release database (db-2026.10.02)."""
    ids = sorted(set(pd.read_csv(os.path.join(BASE_DIR, "output", "pLIN_reference_assignments.tsv"),
                                 sep="\t", usecols=["plasmid_id"]).plasmid_id))
    training = fasta_paths()
    return {i: training.get(i, os.path.join(BASE_DIR, "reference", f"{i}.fasta")) for i in ids}


def _pyrodigal_meta(path, out):
    """Proteins in Prodigal's .faa format, metagenomic mode (pyrodigal, as in the app's plin_backbone.py)."""
    import pyrodigal
    finder = pyrodigal.GeneFinder(meta=True)
    with open(out, "w") as fh:
        for rec in SeqIO.parse(path, "fasta"):
            finder.find_genes(bytes(rec.seq)).write_translations(fh, sequence_id=rec.id)


def _predict_one(args):
    """Returns (plasmid_id, 4-mer vector)."""
    pid, path = args
    faa = os.path.join(FAA, f"{pid}.faa")
    if not os.path.exists(faa):
        _pyrodigal_meta(path, faa + ".part")
        os.replace(faa + ".part", faa)
    return pid, kmer_vector(str(next(SeqIO.parse(path, "fasta")).seq))


def cmd_predict(args):
    os.makedirs(FAA, exist_ok=True)
    plasmids = release_plasmids()
    ids = sorted(plasmids)
    vec = {}
    with ProcessPoolExecutor(args.threads) as ex:
        for n, (pid, v) in enumerate(ex.map(_predict_one, [(i, plasmids[i]) for i in ids], chunksize=16), 1):
            vec[pid] = v
            if n % 5000 == 0:
                print(f"  {n:,}/{len(ids):,}", flush=True)
    np.savez_compressed(os.path.join(OUT, "kmer_vectors.npz"), ids=np.array(ids),
                        X=np.array([vec[i] for i in ids], dtype=np.float32))
    print(f"predicted proteins and 4-mer vectors for {len(ids):,} plasmids")


# ── family catalogue ─────────────────────────────────────────────────────────

def cmd_catalogue(args):
    ids = sorted(release_plasmids())
    combined = os.path.join(OUT, "all_proteins.faa")
    if not os.path.exists(combined):
        n = 0
        with open(combined + ".part", "w") as out:
            for pid in ids:
                for line in open(os.path.join(FAA, f"{pid}.faa")):
                    if line.startswith(">"):
                        out.write(f">{pid}|{line[1:].split()[0]}\n")
                        n += 1
                    else:
                        out.write(line)
        os.replace(combined + ".part", combined)
        print(f"{n:,} proteins from {len(ids):,} plasmids", flush=True)
    prefix = os.path.join(OUT, "clu")
    if not os.path.exists(prefix + "_cluster.tsv"):
        subprocess.run([MMSEQS, "easy-cluster", combined, prefix, os.path.join(OUT, "mmseqs_tmp"),
                        "--min-seq-id", "0.5", "-c", "0.8", "--cov-mode", "0", "--threads", str(args.threads)],
                       check=True)
    clu = pd.read_csv(prefix + "_cluster.tsv", sep="\t", header=None, names=["rep", "member"])
    clu["plasmid_id"] = clu.member.str.split("|").str[0]
    # frozen family IDs: numbered by first appearance in accession order
    order = {pid: k for k, pid in enumerate(ids)}
    clu["rank"] = clu.plasmid_id.map(order)
    first = clu.groupby("rep")["rank"].min().sort_values(kind="stable")
    fam_id = pd.Series(np.arange(1, len(first) + 1), index=first.index)
    clu["family"] = clu.rep.map(fam_id)
    clu[["member", "plasmid_id", "family"]].to_csv(os.path.join(OUT, "protein_families.tsv.gz"), sep="\t",
                                                   index=False, compression="gzip")
    sets = clu.groupby("plasmid_id")["family"].apply(lambda s: np.unique(s.values))
    lens = np.array([len(sets.get(i, [])) for i in ids])
    flat = np.concatenate([sets.get(i, np.array([], dtype=np.int64)) for i in ids]).astype(np.int64)
    np.savez_compressed(os.path.join(OUT, "family_sets.npz"), ids=np.array(ids),
                        offsets=np.concatenate([[0], np.cumsum(lens)]), families=flat,
                        family_rep=fam_id.index.values.astype(str))
    print(f"{len(fam_id):,} families; plasmids without proteins: {(lens == 0).sum():,}")


def cmd_canonicalize(args):
    """Put every identical protein sequence in one family (the lowest ID).

    MMseqs2 clustering occasionally splits identical sequences across families
    (10,309 of 4.47 M unique sequences in db-2026.10.02). Then a protein's
    family depends on which plasmid it came from, and a query plasmid cannot
    reproduce a database plasmid's family set from sequence alone. Deviation 3
    in the pre-registration. The uncanonicalised sets are kept in
    family_sets_raw.npz.
    """
    raw = os.path.join(OUT, "family_sets_raw.npz")
    if not os.path.exists(raw):
        os.replace(os.path.join(OUT, "family_sets.npz"), raw)
    z = np.load(raw, allow_pickle=True)
    ids, reps = z["ids"].tolist(), z["family_rep"]
    fam = pd.read_csv(os.path.join(OUT, "protein_families.tsv.gz"), sep="\t")
    member_hash = {}
    name, seq = None, []
    for line in open(os.path.join(OUT, "all_proteins.faa")):
        if line.startswith(">"):
            if name:
                member_hash[name] = seq_hash("".join(seq))
            name, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    member_hash[name] = seq_hash("".join(seq))
    fam["hash"] = fam.member.map(member_hash).astype(np.uint64)
    fam["canonical"] = fam.groupby("hash")["family"].transform("min")
    changed = int((fam.canonical != fam.family).sum())
    sets = fam.groupby("plasmid_id")["canonical"].apply(lambda s: np.unique(s.values))
    lens = np.array([len(sets.get(i, [])) for i in ids])
    flat = np.concatenate([sets.get(i, np.array([], dtype=np.int64)) for i in ids]).astype(np.int64)
    rep_family = fam.set_index("member").canonical.reindex(reps).values.astype(np.int64)
    np.savez_compressed(os.path.join(OUT, "family_sets.npz"), ids=np.array(ids),
                        offsets=np.concatenate([[0], np.cumsum(lens)]), families=flat,
                        family_rep=reps, family_rep_family=rep_family)
    look = fam.drop_duplicates("hash").sort_values("hash")
    np.savez(EXACT_PROTEINS, hash=look.hash.values, family=look.canonical.values.astype(np.int64))
    print(f"{changed:,} protein copies moved to their sequence's lowest family; "
          f"{fam.family.nunique():,} -> {fam.canonical.nunique():,} families in use")


def load_family_sets():
    z = np.load(os.path.join(OUT, "family_sets.npz"), allow_pickle=True)
    ids, off, fam = z["ids"].tolist(), z["offsets"], z["families"]
    return {i: fam[off[k]:off[k + 1]] for k, i in enumerate(ids)}


def load_vectors():
    z = np.load(os.path.join(OUT, "kmer_vectors.npz"), allow_pickle=True)
    return dict(zip(z["ids"].tolist(), z["X"].astype(np.float64)))


# ── founder tree ─────────────────────────────────────────────────────────────

class _ProteinNode:
    """Founders of one cluster's children, with an inverted index family -> founder positions."""
    __slots__ = ("ids", "sizes", "post")

    def __init__(self):
        self.ids, self.sizes, self.post = [], [], {}

    def first_hit(self, fams, t):
        """Earliest founder with containment >= t, as (position, containment), or (None, None)."""
        if not self.ids or not len(fams):
            return None, None
        hits = [self.post[f] for f in fams.tolist() if f in self.post]
        if not hits:
            return None, None
        inter = np.bincount(np.concatenate(hits), minlength=len(self.ids))
        cont = inter / np.minimum(np.asarray(self.sizes), len(fams))
        ok = np.flatnonzero(cont >= t - 1e-12)
        return (int(ok[0]), float(cont[ok[0]])) if len(ok) else (None, None)

    def append(self, cid, fams):
        j = len(self.ids)
        self.ids.append(cid)
        self.sizes.append(len(fams))
        for f in fams.tolist():
            self.post.setdefault(f, []).append(j)


class V4Tree:
    """Hierarchical founder index: protein containment at L1–L4, 4-mer cosine at L5–L6."""

    def __init__(self, l3):
        self.levels = ["A", "B", "C", "D", "E", "F"]
        self.kinds = ["protein"] * 4 + ["kmer"] * 2
        t = dict(PROTEIN_THRESHOLDS, C=l3)
        self.thresholds = [t["A"], t["B"], t["C"], t["D"], KMER_THRESHOLDS["E"], KMER_THRESHOLDS["F"]]
        self.children = {}
        self.next_id = [1] * 6

    def assign(self, fams, vector, add=True):
        """Return the code tuple; with add=False the tree is unchanged (new clusters get provisional IDs)."""
        u = _unit(vector)
        prefix, provisional = (), list(self.next_id)
        no_protein = len(fams) == 0
        for li, (kind, t) in enumerate(zip(self.kinds, self.thresholds)):
            if kind == "protein" and no_protein:
                prefix += (0,)                             # no proteins: L1–L4 bucket 0
                continue
            node = self.children.get(prefix)
            chosen = None
            if node is not None:
                if kind == "protein":
                    j, _ = node.first_hit(fams, t)
                    chosen = node.ids[j] if j is not None else None
                elif node.n:
                    d = 1.0 - (node.vecs() * u).sum(axis=1)
                    hit = np.flatnonzero(d <= t)
                    chosen = node.ids[int(hit[0])] if len(hit) else None
            if chosen is None:
                if add:
                    chosen = self.next_id[li]
                    self.next_id[li] += 1
                    if node is None:
                        node = self.children[prefix] = _ProteinNode() if kind == "protein" else _Node()
                    node.append(chosen, fams if kind == "protein" else u)
                else:
                    chosen = provisional[li]
                    provisional[li] += 1
            prefix += (chosen,)
        return prefix


def build_codes(ids, fams, vecs, l3):
    """Founder codes for ids in the given (canonical) order; returns (codes dict, tree)."""
    tree = V4Tree(l3)
    codes = {}
    for n, i in enumerate(ids, 1):
        codes[i] = tree.assign(fams[i], vecs[i])
        if n % 20000 == 0:
            print(f"  {n:,}/{len(ids):,}", flush=True)
    return codes, tree


def cmd_build(args):
    if args.l3 not in L3_GRID:
        sys.exit(f"--l3 must be one of the pre-registered grid values {L3_GRID}")
    fams, vecs = load_family_sets(), load_vectors()
    ids = sorted(fams)                                     # canonical order: accession
    codes, _ = build_codes(ids, fams, vecs, args.l3)
    out = os.path.join(OUT, "codes", f"plin_v4_codes_L3_{args.l3:.2f}.tsv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    pd.DataFrame({"plasmid_id": ids, "pLIN_v4": [".".join(map(str, codes[i])) for i in ids],
                  "n_families": [len(fams[i]) for i in ids]}).to_csv(out, sep="\t", index=False)
    print(f"wrote {out}")


# ── release index and query mode ─────────────────────────────────────────────

RELEASE_L3 = 0.40          # calibrated on the calibration half (evaluation/calibration.tsv)
REP_DB = os.path.join(OUT, "mmseqs_reps", "reps")
TREE_FILE = os.path.join(OUT, f"plin_v4_tree_L3_{RELEASE_L3:.2f}.npz")
SEARCH_ARGS = ["--min-seq-id", "0.5", "-c", "0.8", "--cov-mode", "0"]   # same thresholds as the catalogue
EXACT_PROTEINS = os.path.join(OUT, "exact_lookup.npz")       # protein-sequence hash -> family
EXACT_PLASMIDS = os.path.join(OUT, "plasmid_hashes.tsv")     # plasmid-sequence hash -> plasmid_id


def seq_hash(seq):
    """64-bit BLAKE2b of an upper-case sequence (protein: trailing stop removed)."""
    import hashlib
    return int.from_bytes(hashlib.blake2b(seq.upper().rstrip("*").encode(), digest_size=8).digest(), "little")


def _plasmid_hash(args):
    pid, path = args
    return pid, seq_hash("".join(str(r.seq) for r in SeqIO.parse(path, "fasta")))


def save_tree(path, ids, codes, fams, vecs, l3):
    """Compact release tree: one row per cluster with its founder's family set and raw 4-mer vector.

    The founder of a cluster is the first plasmid, in processing order, that carries its code prefix.
    Loading replays the founders in creation order (ascending cluster ID within each node), which repeats
    the exact arithmetic of the build, so every code re-queries identically.
    """
    first = {}
    for i in ids:
        c = codes[i]
        for k in range(6):
            if k < 4 and c[k] == 0:                    # no-protein bucket: not a protein cluster
                continue
            first.setdefault(c[:k + 1], i)
    rows = sorted(first.items(), key=lambda kv: (len(kv[0]), kv[0][-1]))
    level = np.array([len(p) - 1 for p, _ in rows], dtype=np.int8)
    parent = np.array([".".join(map(str, p[:-1])) for p, _ in rows])
    cid = np.array([p[-1] for p, _ in rows], dtype=np.int64)
    founder = np.array([f for _, f in rows])
    prot = [f for (p, f) in rows if len(p) <= 4]
    lens = np.array([len(fams[f]) for f in prot], dtype=np.int64)
    np.savez_compressed(path, level=level, parent=parent, cid=cid, founder=founder,
                        fam_offsets=np.concatenate([[0], np.cumsum(lens)]),
                        fam_values=np.concatenate([fams[f] for f in prot]).astype(np.int64),
                        kmer=np.array([vecs[f] for (p, f) in rows if len(p) > 4], dtype=np.float32),
                        next_id=np.array(V4TreeNext(codes)), l3=np.array(l3))


def V4TreeNext(codes):
    """Next free cluster ID per level after the release (IDs are global per level)."""
    return [max(c[k] for c in codes.values()) + 1 for k in range(6)]


def load_tree(path):
    z = np.load(path, allow_pickle=False)
    tree = V4Tree(float(z["l3"]))
    off, vals, kmer = z["fam_offsets"], z["fam_values"], z["kmer"]
    p_i = k_i = 0
    for lv, par, c in zip(z["level"], z["parent"], z["cid"]):
        prefix = tuple(int(x) for x in str(par).split(".")) if par else ()
        node = tree.children.get(prefix)
        if lv < 4:
            if node is None:
                node = tree.children[prefix] = _ProteinNode()
            node.append(int(c), vals[off[p_i]:off[p_i + 1]])
            p_i += 1
        else:
            if node is None:
                node = tree.children[prefix] = _Node()
            node.append(int(c), _unit(kmer[k_i].astype(np.float64)))
            k_i += 1
    tree.next_id = z["next_id"].tolist()
    return tree


def cmd_index(args):
    """Searchable family-representative database + the compact release founder tree."""
    if not os.path.exists(REP_DB + ".index"):
        os.makedirs(os.path.dirname(REP_DB), exist_ok=True)
        subprocess.run([MMSEQS, "createdb", os.path.join(OUT, "clu_rep_seq.fasta"), REP_DB], check=True)
        subprocess.run([MMSEQS, "createindex", REP_DB, os.path.join(OUT, "mmseqs_reps", "tmp"),
                        "--threads", str(args.threads)], check=True)
    if not os.path.exists(EXACT_PLASMIDS):
        plasmids = release_plasmids()
        with ProcessPoolExecutor(args.threads) as ex:
            hashes = list(ex.map(_plasmid_hash, sorted(plasmids.items()), chunksize=64))
        pd.DataFrame(hashes, columns=["plasmid_id", "hash"]).astype({"hash": str}).to_csv(
            EXACT_PLASMIDS, sep="\t", index=False)
    fams, vecs = load_family_sets(), load_vectors()
    ids = sorted(fams)
    codes, _ = build_codes(ids, fams, vecs, RELEASE_L3)
    save_tree(TREE_FILE, ids, codes, fams, vecs, RELEASE_L3)
    tree = load_tree(TREE_FILE)                        # verify: every plasmid re-queries to its code
    bad = [i for i in ids if tree.assign(fams[i], vecs[i], add=False) != codes[i]]
    if bad:
        sys.exit(f"compact tree does not reproduce {len(bad)} codes, e.g. {bad[:3]}")
    print(f"index at {REP_DB}; compact tree ({os.path.getsize(TREE_FILE) / 1e6:.0f} MB) reproduces all "
          f"{len(ids):,} codes: {TREE_FILE}")


def load_release():
    """(tree, representative name -> family ID, next free family ID)."""
    tree = load_tree(TREE_FILE)
    z = np.load(os.path.join(OUT, "family_sets.npz"), allow_pickle=True)
    reps = z["family_rep"].tolist()
    if "family_rep_family" in z:                       # canonicalised catalogue (deviation 3)
        rep_to_fam = dict(zip(reps, z["family_rep_family"].tolist()))
    else:
        rep_to_fam = {r: k for k, r in enumerate(reps, 1)}
    return tree, rep_to_fam, len(reps) + 1


def query_families(queries, rep_to_fam, next_family, workdir, threads, protein_exact=True):
    """queries: {query_id: fasta path} -> {query_id: sorted family-ID array}.

    Each protein takes the family of its best-scoring representative hit that
    meets the catalogue thresholds; proteins with no hit are clustered among
    themselves at the same thresholds into provisional new families
    (IDs from next_family upwards, not saved to the release).
    """
    os.makedirs(workdir, exist_ok=True)
    faa = os.path.join(workdir, "query.faa")
    owner = {}
    with open(faa, "w") as out:
        for qid, path in queries.items():
            part = os.path.join(workdir, "one.faa")
            _pyrodigal_meta(path, part)
            for line in open(part):
                if line.startswith(">"):
                    tag = f"{qid}|{line[1:].split()[0]}"
                    owner[tag] = qid
                    out.write(f">{tag}\n")
                else:
                    out.write(line)
    fam_of = {}
    z = np.load(EXACT_PROTEINS)
    ex_h, ex_f = z["hash"], z["family"]
    seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(faa, "fasta")}
    for tag, sq in (seqs.items() if protein_exact else ()):   # 1. exact sequence already in the catalogue
        h = np.uint64(seq_hash(sq))
        k = np.searchsorted(ex_h, h)
        if k < len(ex_h) and ex_h[k] == h:
            fam_of[tag] = int(ex_f[k])
    rest = [t for t in owner if t not in fam_of]
    if rest:                                          # 2. MMseqs2 search for unseen proteins
        rfaa = os.path.join(workdir, "unseen.faa")
        with open(rfaa, "w") as out:
            for t in rest:
                out.write(f">{t}\n{seqs[t]}\n")
        faa = rfaa
        m8 = os.path.join(workdir, "hits.m8")
        subprocess.run([MMSEQS, "easy-search", faa, REP_DB, m8, os.path.join(workdir, "tmp"), *SEARCH_ARGS,
                        "--format-output", "query,target,bits", "--threads", str(threads)],
                       check=True, capture_output=True)
        hits = pd.read_csv(m8, sep="\t", header=None, names=["q", "t", "bits"])
        best = hits.sort_values(["q", "bits"], ascending=[True, False], kind="stable").drop_duplicates("q")
        fam_of.update({q: rep_to_fam[t] for q, t in zip(best.q, best.t)})
        novel = [t for t in owner if t not in fam_of]
        if novel:
            nfaa = os.path.join(workdir, "novel.faa")
            keep = set(novel)
            with open(nfaa, "w") as out:
                for rec in SeqIO.parse(faa, "fasta"):
                    if rec.id in keep:
                        out.write(f">{rec.id}\n{rec.seq}\n")
            pre = os.path.join(workdir, "novel")
            subprocess.run([MMSEQS, "easy-cluster", nfaa, pre, os.path.join(workdir, "tmp2"), *SEARCH_ARGS,
                            "--threads", str(threads)], check=True, capture_output=True)
            clu = pd.read_csv(pre + "_cluster.tsv", sep="\t", header=None, names=["rep", "member"])
            new_id = {r: next_family + k for k, r in enumerate(sorted(clu.rep.unique()))}
            fam_of.update({m: new_id[r] for r, m in zip(clu.rep, clu.member)})
    out = {qid: [] for qid in queries}
    for tag, qid in owner.items():
        out[qid].append(fam_of[tag])
    return {q: np.unique(np.array(v, dtype=np.int64)) for q, v in out.items()}


def cmd_query(args):
    """Type FASTA files with the frozen release (same engine as the app: plin_v4_typer)."""
    import time
    from plin_v4_typer import PlinV4Release
    t0 = time.perf_counter()
    rel = PlinV4Release(args.release)
    t_load = time.perf_counter() - t0
    records = []
    for p in (l.strip() for l in open(args.list) if l.strip()):
        recs = list(SeqIO.parse(p, "fasta"))
        records.append((os.path.basename(p).rsplit(".", 1)[0], "".join(str(r.seq) for r in recs)))
    t1 = time.perf_counter()
    df = rel.type_sequences(records, threads=args.threads, workdir=args.workdir,
                            whole_match=not args.no_exact, protein_exact=not args.search_only)
    t_query = time.perf_counter() - t1
    df.rename(columns={"matched_accession": "matched"}).to_csv(args.out, sep="\t", index=False)
    print(f"{len(df)} plasmids typed in {t_query:.1f} s (+ {t_load:.1f} s to load the release); wrote {args.out}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("predict"); p.add_argument("--threads", type=int, default=8)
    p = sub.add_parser("catalogue"); p.add_argument("--threads", type=int, default=8)
    p = sub.add_parser("canonicalize")
    p = sub.add_parser("build"); p.add_argument("--l3", type=float, required=True)
    p = sub.add_parser("index"); p.add_argument("--threads", type=int, default=8)
    p = sub.add_parser("query")
    p.add_argument("list", help="file with one FASTA path per line")
    p.add_argument("--out", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--release", default=os.path.join(OUT, "release"), help="frozen release directory")
    p.add_argument("--no-exact", action="store_true", help="testing: skip whole-plasmid matching")
    p.add_argument("--search-only", action="store_true",
                   help="testing: also skip exact protein matching (worst case: every protein searched)")
    args = ap.parse_args()
    {"predict": cmd_predict, "catalogue": cmd_catalogue, "canonicalize": cmd_canonicalize, "build": cmd_build,
     "index": cmd_index, "query": cmd_query}[args.cmd](args)


if __name__ == "__main__":
    main()
