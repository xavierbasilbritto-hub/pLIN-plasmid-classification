# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
pLIN v4 typing from a frozen release directory (shared by the CLI and the app).

A plasmid gets its code in three steps:
  1. identical to a release plasmid (whole-sequence hash)  -> its published code
  2. otherwise proteins are predicted (pyrodigal, metagenomic mode) and each is
     given a protein family: exact sequence in the catalogue first, then the
     best MMseqs2 hit to a family representative (>= 50% identity, >= 80%
     coverage); proteins with no hit are clustered among themselves into
     provisional families
  3. the founder tree places it: L1–L4 by protein-family containment, L5–L6 by
     4-mer cosine distance within the L4 cluster.
Plasmids typed together extend an in-memory copy of the tree, so they are
coded consistently with each other; any level not in the release is reported
as provisional (column `provisional_from`). The release files are never
modified.

Release directory contents: see build_v4_release.py.
"""

import gzip
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd

from plin_v4 import _pyrodigal_meta, load_tree, seq_hash
from plin_founder import kmer_vector

SEARCH_ARGS = ["--min-seq-id", "0.5", "-c", "0.8", "--cov-mode", "0"]
REQUIRED = ["plin_v4_codes.tsv.gz", "plin_v4_tree.npz", "family_representatives.faa.gz",
            "family_rep_to_family.tsv.gz", "family_exact_lookup.npz", "plasmid_hashes.tsv.gz",
            "DATABASE_VERSION.json"]


def find_mmseqs():
    """Path to an mmseqs binary (PATH, then common conda locations), or None."""
    hit = shutil.which("mmseqs")
    if hit:
        return hit
    for env in ("pLIN_tools", "pLIN_analysis", "base"):
        p = os.path.expanduser(f"~/miniconda3/envs/{env}/bin/mmseqs") if env != "base" \
            else os.path.expanduser("~/miniconda3/bin/mmseqs")
        if os.path.exists(p):
            return p
    return None


def missing_files(release_dir):
    return [f for f in REQUIRED if not os.path.exists(os.path.join(release_dir, f))]


class PlinV4Release:
    def __init__(self, release_dir, mmseqs=None, index_dir=None):
        miss = missing_files(release_dir)
        if miss:
            raise FileNotFoundError(f"pLIN v4 release incomplete in {release_dir}: missing {miss}")
        self.dir = release_dir
        self.mmseqs = mmseqs or find_mmseqs()
        self.index_dir = index_dir or os.path.join(release_dir, "mmseqs_index")
        self.tree_path = os.path.join(release_dir, "plin_v4_tree.npz")
        self.release_next = list(load_tree(self.tree_path).next_id)
        codes = pd.read_csv(os.path.join(release_dir, "plin_v4_codes.tsv.gz"), sep="\t", dtype=str)
        self.codes = codes.set_index("accession")
        h = pd.read_csv(os.path.join(release_dir, "plasmid_hashes.tsv.gz"), sep="\t", dtype=str)
        h["accession"] = h.plasmid_id.str.replace("^RefSeq_", "", regex=True)
        self.known = dict(zip(h.hash.astype("uint64"), h.accession))
        z = np.load(os.path.join(release_dir, "family_exact_lookup.npz"))
        self.ex_hash, self.ex_fam = z["hash"], z["family"]
        r = pd.read_csv(os.path.join(release_dir, "family_rep_to_family.tsv.gz"), sep="\t")
        self.rep_to_fam = dict(zip(r.representative, r.family))
        self.next_family = int(r.family.max()) + 1

    # ── MMseqs2 search target ────────────────────────────────────────────────
    def search_db(self, threads=4):
        """Build (once) the MMseqs2 search target; returns (path, kind).

        Preferred target: every unique catalogued protein (family_members.faa.gz, headers = sequence
        hash), so a protein takes the family of its nearest catalogued protein; a point mutation
        then keeps its original's family. Fallback: one representative per family.
        """
        members = os.path.join(self.dir, "family_members.faa.gz")
        kind = "members" if os.path.exists(members) else "representatives"
        db = os.path.join(self.index_dir, kind)
        if os.path.exists(db + ".done"):
            return db, kind
        if not self.mmseqs:
            raise RuntimeError("MMseqs2 not found: install it (conda install -c bioconda mmseqs2) for pLIN v4")
        os.makedirs(self.index_dir, exist_ok=True)
        faa = os.path.join(self.index_dir, f"{kind}.faa")
        src_f = members if kind == "members" else os.path.join(self.dir, "family_representatives.faa.gz")
        with gzip.open(src_f, "rb") as src, open(faa, "wb") as dst:
            shutil.copyfileobj(src, dst)
        subprocess.run([self.mmseqs, "createdb", faa, db], check=True, capture_output=True)
        subprocess.run([self.mmseqs, "createindex", db, os.path.join(self.index_dir, "tmp"), "--threads", str(threads)],
                       check=True, capture_output=True)
        os.remove(faa)
        open(db + ".done", "w").close()
        return db, kind

    # ── typing ───────────────────────────────────────────────────────────────
    def _families(self, fasta_by_id, workdir, threads, protein_exact=True):
        faa = os.path.join(workdir, "query.faa")
        owner, seqs = {}, {}
        with open(faa, "w") as out:
            for qid, path in fasta_by_id.items():
                part = os.path.join(workdir, "one.faa")
                _pyrodigal_meta(path, part)
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
        if protein_exact:
            for tag, sq in seqs.items():
                h = np.uint64(seq_hash(sq))
                k = np.searchsorted(self.ex_hash, h)
                if k < len(self.ex_hash) and self.ex_hash[k] == h:
                    fam_of[tag] = int(self.ex_fam[k])
        rest = [t for t in owner if t not in fam_of]
        if rest:
            db, kind = self.search_db(threads)
            rfaa = os.path.join(workdir, "unseen.faa")
            with open(rfaa, "w") as out:
                for t in rest:
                    out.write(f">{t}\n{seqs[t]}\n")
            m8 = os.path.join(workdir, "hits.m8")
            subprocess.run([self.mmseqs, "easy-search", rfaa, db, m8, os.path.join(workdir, "tmp"), *SEARCH_ARGS,
                            "--format-output", "query,target,bits", "--threads", str(threads)],
                           check=True, capture_output=True)
            if os.path.getsize(m8):
                hits = pd.read_csv(m8, sep="\t", header=None, names=["q", "t", "bits"], dtype={"q": str, "t": str})
                best = hits.sort_values(["q", "bits"], ascending=[True, False], kind="stable").drop_duplicates("q")
                if kind == "members":
                    hh = np.array([int(x) for x in best.t], dtype=np.uint64)   # exact: no float parsing
                    k = np.searchsorted(self.ex_hash, hh)
                    assert (self.ex_hash[np.minimum(k, len(self.ex_hash) - 1)] == hh).all(), "unknown member hash"
                    fam_of.update({q: int(self.ex_fam[j]) for q, j in zip(best.q, k)})
                else:
                    fam_of.update({q: int(self.rep_to_fam[t]) for q, t in zip(best.q, best.t)})
            novel = [t for t in rest if t not in fam_of]
            if novel:
                nfaa = os.path.join(workdir, "novel.faa")
                with open(nfaa, "w") as out:
                    for t in novel:
                        out.write(f">{t}\n{seqs[t]}\n")
                pre = os.path.join(workdir, "novel")
                subprocess.run([self.mmseqs, "easy-cluster", nfaa, pre, os.path.join(workdir, "tmp2"), *SEARCH_ARGS,
                                "--threads", str(threads)], check=True, capture_output=True)
                clu = pd.read_csv(pre + "_cluster.tsv", sep="\t", header=None, names=["rep", "member"])
                new_id = {r: self.next_family + k for k, r in enumerate(sorted(clu.rep.unique()))}
                fam_of.update({m: new_id[r] for r, m in zip(clu.rep, clu.member)})
        out = {qid: [] for qid in fasta_by_id}
        for tag, qid in owner.items():
            out[qid].append(fam_of[tag])
        return {q: np.unique(np.array(v, dtype=np.int64)) for q, v in out.items()}

    def type_sequences(self, records, threads=4, workdir=None, whole_match=True, protein_exact=True):
        """records: list of (plasmid_id, sequence). Returns a DataFrame, one row per record, input order."""
        tmp = workdir or tempfile.mkdtemp(prefix="plin_v4_")
        os.makedirs(tmp, exist_ok=True)
        rows, todo = {}, {}
        for pid, seq in records:
            match = self.known.get(np.uint64(seq_hash(seq))) if whole_match else None
            if match is not None and match in self.codes.index:
                r = self.codes.loc[match]
                rows[pid] = {"pLIN_v4": r.pLIN_v4, "matched_accession": match, "provisional_from": None,
                             "n_families": int(r.n_families), "pLIN_v3_of_match": r.pLIN_v3}
            else:
                path = os.path.join(tmp, f"q{len(todo)}.fasta")
                with open(path, "w") as fh:
                    fh.write(f">{pid}\n{seq}\n")
                todo[pid] = path
        if todo:
            fams = self._families(todo, tmp, threads, protein_exact=protein_exact)
            tree = load_tree(self.tree_path)  # fresh copy per call: earlier calls never leak in
            for pid in sorted(todo):          # deterministic order within a session
                seq = next(s for i, s in records if i == pid)
                code = tree.assign(fams[pid], kmer_vector(seq), add=True)
                new = next((f"L{k + 1}" for k, c in enumerate(code) if c >= self.release_next[k]), None)
                rows[pid] = {"pLIN_v4": ".".join(map(str, code)), "matched_accession": None,
                             "provisional_from": new, "n_families": int(len(fams[pid])), "pLIN_v3_of_match": None}
        if workdir is None:
            shutil.rmtree(tmp, ignore_errors=True)
        df = pd.DataFrame([{"plasmid_id": pid, **rows[pid]} for pid, _ in records])
        lv = df.pLIN_v4.str.split(".", expand=True)
        for k in range(6):
            df[f"L{k + 1}"] = lv[k].astype(int)
        return df
