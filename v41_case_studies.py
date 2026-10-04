#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
pLIN v4.1 on the two DEVELOPMENT case studies (used during design; descriptive only).

  Swiss VIM-1  19 plasmid contigs from 8 isolates, typed in one session with
               plin_v41_typer; for the 8 blaVIM-1 plasmids, levels shared per pair
               vs the alignment (output/backbone_v4/swiss_vim1/vim_pairs.tsv)
  outbreak     74 published outbreak plasmids (27 studies), typed in one session;
               same-study and different-study pairs grouped per level

Usage:
  python v41_case_studies.py [--index-dir DIR]
Output: output/backbone_v41/case_studies/{swiss_vim1_codes.tsv, swiss_vim1_pairs.tsv,
        outbreak_codes.tsv, outbreak_metrics.tsv}
"""

import argparse
import glob
import os
from itertools import combinations

import numpy as np
import pandas as pd
from Bio import SeqIO

from plin_v41_typer import PlinV41Release

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "output", "backbone_v41", "case_studies")


def shared(a, b):
    n = 0
    for x, y in zip(a.split("."), b.split(".")):
        if x != y:
            break
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index-dir", default=None)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rel = PlinV41Release(os.path.join(BASE_DIR, "output", "backbone_v41", "release"), index_dir=args.index_dir)

    recs = []
    for f in sorted(glob.glob(os.path.join(BASE_DIR, "sample_data", "swiss_vim1_outbreak", "*_plasmids.fasta"))):
        recs += [(r.id, str(r.seq)) for r in SeqIO.parse(f, "fasta")]
    vim = rel.type_sequences(recs, threads=8)
    prev = pd.read_csv(os.path.join(BASE_DIR, "output", "backbone_v4", "swiss_vim1", "contig_codes.tsv"), sep="\t") \
        .set_index("contig_id")
    for col in ("pLIN_single_linkage_v3.2.3", "pLIN_v3_founder", "pLIN_v4", "key_genes", "length_bp", "topology"):
        vim[col] = vim.plasmid_id.map(prev[col])
    vim.to_csv(os.path.join(OUT, "swiss_vim1_codes.tsv"), sep="\t", index=False)
    code = dict(zip(vim.plasmid_id, vim.pLIN_v41))
    vp = pd.read_csv(os.path.join(BASE_DIR, "output", "backbone_v4", "swiss_vim1", "vim_pairs.tsv"), sep="\t")
    vp["levels_shared_v41"] = [shared(code[a], code[b]) for a, b in zip(vp.contig_shorter, vp.contig_longer)]
    vp.to_csv(os.path.join(OUT, "swiss_vim1_pairs.tsv"), sep="\t", index=False)
    s, o = vp[vp.same_lineage_by_alignment], vp[~vp.same_lineage_by_alignment]
    print(f"Swiss VIM-1: same-lineage pairs sharing L5: {(s.levels_shared_v41 >= 5).sum()}/{len(s)}; "
          f"other pairs split before L5: {(o.levels_shared_v41 < 5).sum()}/{len(o)}; "
          f"all blaVIM-1 pairs share L1: {(vp.levels_shared_v41 >= 1).sum()}/{len(vp)}")
    print(vim[["plasmid_id", "length_bp", "key_genes", "pLIN_v41", "provisional_from", "nearest_database_plasmid",
               "nearest_similarity"]].to_string(index=False))

    ob = pd.read_csv(os.path.join(BASE_DIR, "output", "outbreak_validation_founder_results.tsv"), sep="\t")
    files = {os.path.basename(f).rsplit(".", 1)[0]: f for f in
             glob.glob(os.path.join(BASE_DIR, "outbreak_validation", "*.fasta")) +
             glob.glob(os.path.join(BASE_DIR, "outbreak_validation", "expanded_sequences", "*.fasta"))}
    recs = [(a, "".join(str(r.seq) for r in SeqIO.parse(files[a], "fasta"))) for a in ob.accession]
    oc = rel.type_sequences(recs, threads=8)
    oc["study"] = oc.plasmid_id.map(dict(zip(ob.accession, ob.study)))
    oc.to_csv(os.path.join(OUT, "outbreak_codes.tsv"), sep="\t", index=False)
    c, st = dict(zip(oc.plasmid_id, oc.pLIN_v41)), dict(zip(oc.plasmid_id, oc.study))
    rows = []
    for k in range(1, 7):
        same, diff = [], []
        for a, b in combinations(oc.plasmid_id, 2):
            (same if st[a] == st[b] else diff).append(shared(c[a], c[b]) >= k)
        rows.append({"level": f"L{k}", "same_study_together_pct": round(100 * np.mean(same), 1),
                     "different_study_together_pct": round(100 * np.mean(diff), 2), "false_links": int(np.sum(diff))})
    om = pd.DataFrame(rows)
    om.to_csv(os.path.join(OUT, "outbreak_metrics.tsv"), sep="\t", index=False)
    print("\nOutbreak set (74 plasmids, 27 studies):\n" + om.to_string(index=False))


if __name__ == "__main__":
    main()
