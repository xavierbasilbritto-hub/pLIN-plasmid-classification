# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Backbone protein-family comparison of plasmids (companion to pLIN codes).

pLIN codes come from tetranucleotide composition, which identifies
near-identical plasmids well but says little about backbone relatedness at the
coarse levels (L1–L4), and misses some same-lineage pairs where one plasmid
contains the other plus extra DNA. Shared backbone protein families answer
both questions (pilot_protein_backbone.py, 3,183 alignment-validated pairs):

  * related backbone: containment >= 0.5 of the smaller plasmid's backbone
    families detects pairs sharing >= 50% backbone sequence with 92.5%
    sensitivity and 99.0% specificity (AUC 0.980 vs 0.937 for cosine distance;
    0.943 vs 0.817 among pairs sharing only L1–L4);
  * possible same lineage: different L6 code but containment >= 0.9 and cosine
    distance <= 0.005 recovers about half of the same-lineage pairs that L6
    misses (recall 84.5% -> 92.2%, precision 74.3% -> 73.0%).

This is an annotation layer only: it never changes pLIN codes, whose
stability depends on the composition-based founder rule alone.

Proteins are predicted with pyrodigal (metagenomic mode) and clustered into
families with MMseqs2 easy-cluster (>= 50% identity, >= 80% coverage).
Proteins overlapping annotated AMR genes (and any supplied IS/transposon
intervals) by at least half their length are cargo; all others are backbone.
"""

import os
import shutil
import subprocess
import tempfile
from itertools import combinations

import numpy as np
import pandas as pd

RELATED_BACKBONE_MIN = 0.5
SAME_LINEAGE_CONTAINMENT_MIN = 0.9
SAME_LINEAGE_COSINE_MAX = 0.005


def find_mmseqs():
    for cand in (shutil.which("mmseqs"),
                 os.path.expanduser("~/miniconda3/envs/pLIN_tools/bin/mmseqs")):
        if cand and os.path.exists(cand):
            return cand
    return None


def predict_proteins(records):
    """records: list of dicts with plasmid_id and sequence. Returns {pid: [(prot_id, start, end, aa)]}."""
    import pyrodigal
    finder = pyrodigal.GeneFinder(meta=True)
    out = {}
    for rec in records:
        genes = finder.find_genes(rec["sequence"].encode())
        out[rec["plasmid_id"]] = [(f"{rec['plasmid_id']}|{i}", g.begin, g.end, g.translate().rstrip("*"))
                                  for i, g in enumerate(genes, 1)]
    return out


def cargo_intervals_from_amr(amr_df):
    """{plasmid_id: [(start, end)]} from an AMRFinderPlus table (AMR elements only)."""
    if amr_df is None or len(amr_df) == 0:
        return {}
    cols = {c.lower(): c for c in amr_df.columns}
    cid = cols.get("contig id") or cols.get("plasmid_id") or cols.get("source_file")
    s, e = cols.get("start"), cols.get("stop") or cols.get("end")
    typ = cols.get("type") or cols.get("element type")
    if not (cid and s and e):
        return {}
    df = amr_df[amr_df[typ] == "AMR"] if typ else amr_df
    out = {}
    for pid, a, b in zip(df[cid].astype(str), df[s], df[e]):
        out.setdefault(pid, []).append((int(min(a, b)), int(max(a, b))))
    return out


def protein_families(proteins, threads=4, mmseqs=None, min_id=0.5, cov=0.8):
    """Cluster all proteins; returns {prot_id: family_rep}."""
    mmseqs = mmseqs or find_mmseqs()
    if mmseqs is None:
        raise FileNotFoundError("MMseqs2 not found (install with: conda install -c bioconda mmseqs2)")
    with tempfile.TemporaryDirectory() as tmp:
        faa = os.path.join(tmp, "proteins.faa")
        with open(faa, "w") as fh:
            for plist in proteins.values():
                for pid, _, _, aa in plist:
                    if aa:
                        fh.write(f">{pid}\n{aa}\n")
        prefix = os.path.join(tmp, "clu")
        subprocess.run([mmseqs, "easy-cluster", faa, prefix, os.path.join(tmp, "work"),
                        "--min-seq-id", str(min_id), "-c", str(cov), "--threads", str(threads)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        clu = pd.read_csv(prefix + "_cluster.tsv", sep="\t", header=None, names=["rep", "member"])
    return dict(zip(clu.member, clu.rep))


def backbone_sharing(records, amr_df=None, extra_cargo=None, threads=4, mmseqs=None):
    """Pairwise backbone protein-family sharing for a set of plasmids.

    Returns a DataFrame with one row per pair: backbone family counts, shared
    families, containment (shared / families of the plasmid with fewer) and
    Jaccard, plus the related-backbone call.
    """
    proteins = predict_proteins(records)
    cargo = cargo_intervals_from_amr(amr_df)
    for pid, iv in (extra_cargo or {}).items():
        cargo.setdefault(pid, []).extend(iv)
    fam = protein_families(proteins, threads=threads, mmseqs=mmseqs)

    backbone = {}
    for pid, plist in proteins.items():
        ivs = cargo.get(pid, [])
        fams = set()
        for prot, s, e, aa in plist:
            if prot not in fam:
                continue
            is_cargo = any(min(e, b) - max(s, a) + 1 >= 0.5 * (e - s + 1) for a, b in ivs)
            if not is_cargo:
                fams.add(fam[prot])
        backbone[pid] = fams

    rows = []
    for a, b in combinations([r["plasmid_id"] for r in records], 2):
        fa, fb = backbone.get(a, set()), backbone.get(b, set())
        small = min(len(fa), len(fb))
        shared = len(fa & fb)
        containment = shared / small if small else np.nan
        rows.append({"plasmid_1": a, "plasmid_2": b,
                     "backbone_families_1": len(fa), "backbone_families_2": len(fb),
                     "shared_families": shared,
                     "containment": round(containment, 3) if small else np.nan,
                     "jaccard": round(shared / len(fa | fb), 3) if (fa | fb) else np.nan,
                     "related_backbone": bool(small and containment >= RELATED_BACKBONE_MIN)})
    return pd.DataFrame(rows)


def flag_possible_same_lineage(pairs, codes, cosine):
    """Add the 'possible same lineage' flag: different L6 code, but backbone
    containment >= 0.9 and cosine distance <= 0.005.

    codes: {plasmid_id: pLIN code}; cosine: function (id1, id2) -> distance.
    """
    out = pairs.copy()
    out["same_L6"] = [codes.get(a) == codes.get(b) for a, b in zip(out.plasmid_1, out.plasmid_2)]
    out["cosine_distance"] = [round(float(cosine(a, b)), 5) for a, b in zip(out.plasmid_1, out.plasmid_2)]
    out["possible_same_lineage"] = (~out.same_L6) & (out.containment >= SAME_LINEAGE_CONTAINMENT_MIN) & \
                                   (out.cosine_distance <= SAME_LINEAGE_COSINE_MAX)
    return out
