# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
pLIN v4.1 typing from a frozen release directory (used by the CLI and the app).

For each plasmid sequence:
  1. identical to a release plasmid (whole-sequence hash)  -> its published code
  2. otherwise: proteins (pyrodigal, metagenomic mode) get catalogue families:
     identical sequence first, else the family of the nearest catalogued protein
     (MMseqs2, >= 50% identity, >= 80% coverage), else a new provisional family:
     and a k-mer sketch is taken (plin_kmers)
  3. plin_v41_nn.HybridIndex assigns the code: a plasmid with a database (or
     same-session) relative at symmetric k-mer similarity >= 0.80 copies that
     relative's backbone and lineage; otherwise L1–L4 come from the founder tree
     and L5–L6 are new.
Plasmids typed together form one session: they can group with each other, new
clusters exist only in the session (column `provisional_from`), and the release
files are never modified. Each row also reports the most similar database
plasmid and its symmetric k-mer similarity.

Release directory: build_v41_release.py.
"""

import gzip
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd

from plin_kmers import adaptive_sketch, min_containment
from plin_v4 import _pyrodigal_meta, seq_hash
from plin_v41_nn import load_hybrid

SEARCH_ARGS = ["--min-seq-id", "0.5", "-c", "0.8", "--cov-mode", "0"]
REQUIRED = ["plin_v41_codes.tsv.gz", "plin_v41_index.npz", "family_members.faa.gz", "family_exact_lookup.npz",
            "plasmid_hashes.tsv.gz", "DATABASE_VERSION.json"]
LEVEL_NAMES = ["backbone family", "backbone group", "shared backbone", "backbone variant", "lineage",
               "near-identical"]


def find_mmseqs():
    hit = shutil.which("mmseqs")
    if hit:
        return hit
    for p in ("~/miniconda3/envs/pLIN_tools/bin/mmseqs", "~/miniconda3/bin/mmseqs", "~/miniforge3/bin/mmseqs"):
        if os.path.exists(os.path.expanduser(p)):
            return os.path.expanduser(p)
    return None


def missing_files(release_dir):
    return [f for f in REQUIRED if not os.path.exists(os.path.join(release_dir, f))]


class PlinV41Release:
    def __init__(self, release_dir, mmseqs=None, index_dir=None):
        miss = missing_files(release_dir)
        if miss:
            raise FileNotFoundError(f"pLIN v4.1 release incomplete in {release_dir}: missing {miss}")
        self.dir = release_dir
        self.mmseqs = mmseqs or find_mmseqs()
        self.index_dir = index_dir or os.path.join(release_dir, "mmseqs_index")
        self.hx = load_hybrid(os.path.join(release_dir, "plin_v41_index.npz"))
        self.n = len(self.hx.codes)
        codes = pd.read_csv(os.path.join(release_dir, "plin_v41_codes.tsv.gz"), sep="\t", dtype=str)
        self.accessions = codes.accession.tolist()
        self.codes = codes.set_index("accession")
        h = pd.read_csv(os.path.join(release_dir, "plasmid_hashes.tsv.gz"), sep="\t", dtype=str)
        self.known = dict(zip(h.hash.astype("uint64"), h.accession))
        z = np.load(os.path.join(release_dir, "family_exact_lookup.npz"))
        self.ex_hash, self.ex_fam = z["hash"], z["family"]
        self.next_family = int(self.ex_fam.max()) + 1
        self.release_next = list(self.hx.tree.next_id) + list(self.hx.next_lin)

    def search_db(self, threads=4):
        """Build (once) the MMseqs2 database of all catalogued proteins (~16 GB with index)."""
        db = os.path.join(self.index_dir, "members")
        if os.path.exists(db + ".done"):
            return db
        if not self.mmseqs:
            raise RuntimeError("MMseqs2 not found: install it (conda install -c bioconda mmseqs2) to type new plasmids")
        os.makedirs(self.index_dir, exist_ok=True)
        faa = os.path.join(self.index_dir, "members.faa")
        with gzip.open(os.path.join(self.dir, "family_members.faa.gz"), "rb") as src, open(faa, "wb") as dst:
            shutil.copyfileobj(src, dst)
        subprocess.run([self.mmseqs, "createdb", faa, db], check=True, capture_output=True)
        subprocess.run([self.mmseqs, "createindex", db, os.path.join(self.index_dir, "tmp"), "--threads", str(threads)],
                       check=True, capture_output=True)
        os.remove(faa)
        open(db + ".done", "w").close()
        return db

    def families(self, records, workdir, threads=4):
        """records: list of (id, sequence) -> {id: sorted family-ID array}."""
        seqs, owner = {}, {}
        for qid, seq in records:
            fa = os.path.join(workdir, "one.fa")
            with open(fa, "w") as fh:
                fh.write(f">{qid}\n{seq}\n")
            part = os.path.join(workdir, "one.faa")
            _pyrodigal_meta(fa, part)
            tag, buf = None, []
            for line in list(open(part)) + [">"]:
                if line.startswith(">"):
                    if tag:
                        seqs[tag] = "".join(buf)
                    if line.strip() == ">":
                        break
                    tag, buf = f"{qid}|{line[1:].split()[0]}", []
                    owner[tag] = qid
                else:
                    buf.append(line.strip())
        fam_of = {}
        for tag, sq in seqs.items():
            hh = np.uint64(seq_hash(sq))
            k = np.searchsorted(self.ex_hash, hh)
            if k < len(self.ex_hash) and self.ex_hash[k] == hh:
                fam_of[tag] = int(self.ex_fam[k])
        rest = [t for t in owner if t not in fam_of]
        if rest:
            db = self.search_db(threads)
            rfaa, m8 = os.path.join(workdir, "unseen.faa"), os.path.join(workdir, "hits.m8")
            with open(rfaa, "w") as out:
                for t in rest:
                    out.write(f">{t}\n{seqs[t]}\n")
            subprocess.run([self.mmseqs, "easy-search", rfaa, db, m8, os.path.join(workdir, "tmp"), *SEARCH_ARGS,
                            "--format-output", "query,target,bits", "--threads", str(threads)],
                           check=True, capture_output=True)
            if os.path.getsize(m8):
                hits = pd.read_csv(m8, sep="\t", header=None, names=["q", "t", "bits"], dtype={"q": str, "t": str})
                best = hits.sort_values(["q", "bits"], ascending=[True, False], kind="stable").drop_duplicates("q")
                hh = np.array([int(x) for x in best.t], dtype=np.uint64)
                k = np.searchsorted(self.ex_hash, hh)
                fam_of.update({q: int(self.ex_fam[j]) for q, j in zip(best.q, k)})
            novel = [t for t in rest if t not in fam_of]
            if novel:
                nfaa, pre = os.path.join(workdir, "novel.faa"), os.path.join(workdir, "novel")
                with open(nfaa, "w") as out:
                    for t in novel:
                        out.write(f">{t}\n{seqs[t]}\n")
                subprocess.run([self.mmseqs, "easy-cluster", nfaa, pre, os.path.join(workdir, "tmp2"), *SEARCH_ARGS,
                                "--threads", str(threads)], check=True, capture_output=True)
                clu = pd.read_csv(pre + "_cluster.tsv", sep="\t", header=None, names=["rep", "member"])
                new_id = {r: self.next_family + k for k, r in enumerate(sorted(clu.rep.unique()))}
                fam_of.update({m: new_id[r] for r, m in zip(clu.rep, clu.member)})
        out = {qid: [] for qid, _ in records}
        for tag, qid in owner.items():
            out[qid].append(fam_of[tag])
        return {q: np.unique(np.array(v, dtype=np.int64)) for q, v in out.items()}

    def type_sequences(self, records, threads=4, workdir=None, whole_match=True):
        """records: list of (plasmid_id, sequence). Returns a DataFrame in input order."""
        tmp = workdir or tempfile.mkdtemp(prefix="plin_v41_")
        os.makedirs(tmp, exist_ok=True)
        rows, todo = {}, []
        for pid, seq in records:
            seq = seq.upper()
            match = self.known.get(np.uint64(seq_hash(seq))) if whole_match else None
            if match is not None and match in self.codes.index:
                rows[pid] = {"pLIN_v41": self.codes.at[match, "pLIN_v41"], "matched_accession": match,
                             "provisional_from": None, "nearest_database_plasmid": match,
                             "nearest_similarity": 1.0}
            else:
                todo.append((pid, seq))
        if todo:
            fams = self.families(todo, tmp, threads)
            sess = self.hx.session()
            s_sk, s_sc, s_codes = [], [], []
            for pid, seq in sorted(todo):          # deterministic order within a session
                sk, sc = adaptive_sketch(seq)
                extra = None
                if s_codes:
                    extra = {"codes": np.array(s_codes),
                             "kmin": np.array([min_containment(sk, a, sc, b) for a, b in zip(s_sk, s_sc)])}
                code, _ = sess.assign(fams[pid], sk, sc, self.n, extra=extra, add=True)
                rows_, vals = sess.ix.similarities(np.zeros(0, np.int64), sk, sc, self.n)["kmin"]
                j = int(np.argmax(vals)) if len(vals) else -1
                new = next((f"L{k + 1}" for k, c in enumerate(code) if c >= self.release_next[k]), None)
                rows[pid] = {"pLIN_v41": ".".join(map(str, code)), "matched_accession": None, "provisional_from": new,
                             "nearest_database_plasmid": self.accessions[rows_[j]] if j >= 0 else None,
                             "nearest_similarity": round(float(vals[j]), 4) if j >= 0 else 0.0}
                s_sk.append(sk); s_sc.append(sc); s_codes.append(code)
        if workdir is None:
            shutil.rmtree(tmp, ignore_errors=True)
        df = pd.DataFrame([{"plasmid_id": pid, **rows[pid]} for pid, _ in records])
        lv = df.pLIN_v41.str.split(".", expand=True)
        for k in range(6):
            df[f"L{k + 1}"] = lv[k].astype(int)
        return df
