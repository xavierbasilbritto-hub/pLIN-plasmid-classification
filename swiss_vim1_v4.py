#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Swiss VIM-1 outbreak case study across pLIN versions.

The 19 plasmid contigs of the eight ONT-sequenced NARACHVIM isolates
(sample_data/swiss_vim1_outbreak/) are typed with
  single-linkage v3.2.3  codes from the earlier manuscript run
                         (outbreak_validation/swiss_vim1/plin_results/)
  v3 founder             data/plin_founder_tree_reference.npz, session copy
  v4                     frozen release (plin_v4_typer)
AMRFinderPlus finds blaVIM and other key genes, and every pair of
blaVIM-carrying contigs is aligned exactly as in the evaluation
(validate_alignment_backbone.align_pair) to check, independently of any
code, whether they are one lineage (AF_min >= 0.8 and identity >= 99%).

Usage:
  python swiss_vim1_v4.py
Output: output/backbone_v4/swiss_vim1/{contig_codes.tsv, vim_pairs.tsv}
"""

import glob
import os
import shutil
import subprocess
from itertools import combinations

import pandas as pd
from Bio import SeqIO

from plin_founder import FounderTree, kmer_vector, session_copy
from plin_v4_typer import PlinV4Release
from validate_alignment_backbone import align_pair

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SD = os.path.join(BASE_DIR, "sample_data", "swiss_vim1_outbreak")
OUT = os.path.join(BASE_DIR, "output", "backbone_v4", "swiss_vim1")
OLD = os.path.join(BASE_DIR, "outbreak_validation", "swiss_vim1", "plin_results", "swiss_vim1_plin_results.tsv")
KEY = ("blaVIM", "blaNDM", "blaKPC", "blaOXA-48", "blaCTX-M", "mcr", "blaSHV-12", "qnr")


def main():
    os.makedirs(OUT, exist_ok=True)
    recs = []
    for f in sorted(glob.glob(os.path.join(SD, "*_plasmids.fasta"))):
        for r in SeqIO.parse(f, "fasta"):
            topo = r.description.split("topology=")[-1] if "topology=" in r.description else ""
            recs.append({"contig_id": r.id, "sample": r.id.split("_")[0], "length_bp": len(r.seq),
                         "topology": topo, "seq": str(r.seq).upper()})
    recs.sort(key=lambda x: x["contig_id"])
    df = pd.DataFrame(recs)

    old = pd.read_csv(OLD, sep="\t").set_index("contig_id")
    df["pLIN_single_linkage_v3.2.3"] = df.contig_id.map(old.pLIN)

    tree = session_copy(FounderTree.from_npz(os.path.join(BASE_DIR, "data", "plin_founder_tree_reference.npz")))
    df["pLIN_v3_founder"] = [".".join(map(str, tree.assign(kmer_vector(s), key=c)[0])) for c, s in zip(df.contig_id, df.seq)]

    rel = PlinV4Release(os.path.join(BASE_DIR, "output", "backbone_v4", "release"))
    v4 = rel.type_sequences(list(zip(df.contig_id, df.seq)), threads=8).set_index("plasmid_id")
    df["pLIN_v4"] = df.contig_id.map(v4.pLIN_v4)
    df["v4_provisional_from"] = df.contig_id.map(v4.provisional_from)

    # AMRFinderPlus (nucleotide mode)
    fa = os.path.join(OUT, "contigs.fasta")
    with open(fa, "w") as fh:
        for c, s in zip(df.contig_id, df.seq):
            fh.write(f">{c}\n{s}\n")
    amr_bin = shutil.which("amrfinder") or os.path.expanduser("~/miniconda3/envs/pLIN_analysis/bin/amrfinder")
    amr_out = os.path.join(OUT, "amrfinder.tsv")
    subprocess.run([amr_bin, "-n", fa, "--threads", "8", "-o", amr_out], check=True, capture_output=True)
    amr = pd.read_csv(amr_out, sep="\t")
    sym = amr.columns[[c.lower().startswith(("element symbol", "gene symbol")) for c in amr.columns]][0]
    genes = amr.groupby("Contig id")[sym].apply(lambda s: ";".join(sorted(set(s))))
    df["AMR_genes"] = df.contig_id.map(genes).fillna("")
    df["blaVIM"] = df.AMR_genes.str.contains("blaVIM")
    df["key_genes"] = df.AMR_genes.map(lambda g: ";".join(x for x in g.split(";") if x.startswith(KEY)))
    df.drop(columns="seq").to_csv(os.path.join(OUT, "contig_codes.tsv"), sep="\t", index=False)

    # pairwise alignment of blaVIM contigs
    seqs = dict(zip(df.contig_id, df.seq))
    vim = df[df.blaVIM].contig_id.tolist()
    rows = []
    for a, b in combinations(vim, 2):
        a, b = sorted((a, b), key=lambda x: len(seqs[x]))
        pa, pb = os.path.join(OUT, f"{a}.fa"), os.path.join(OUT, f"{b}.fa")
        for p, x in ((pa, a), (pb, b)):
            if not os.path.exists(p):
                open(p, "w").write(f">{x}\n{seqs[x]}\n")
        r = align_pair((pa, pb, len(seqs[a]), len(seqs[b]), [], [], False))
        code = df.set_index("contig_id")
        shared = lambda col: sum(1 for _ in __import__("itertools").takewhile(
            lambda t: t[0] == t[1], zip(str(code.at[a, col]).split("."), str(code.at[b, col]).split("."))))
        rows.append({"contig_shorter": a, "contig_longer": b, "AF_shorter": round(r["a"]["AF"], 3),
                     "AF_longer": round(r["b"]["AF"], 3), "AF_min": round(min(r["a"]["AF"], r["b"]["AF"]), 3),
                     "ANI_aln": round(r["ANI_aln"], 3),
                     "same_lineage_by_alignment": min(r["a"]["AF"], r["b"]["AF"]) >= 0.8 and r["ANI_aln"] >= 99,
                     "levels_shared_single_linkage": shared("pLIN_single_linkage_v3.2.3"),
                     "levels_shared_v3_founder": shared("pLIN_v3_founder"),
                     "levels_shared_v4": shared("pLIN_v4")})
    pairs = pd.DataFrame(rows)
    pairs.to_csv(os.path.join(OUT, "vim_pairs.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["seq", "AMR_genes"]).to_string(index=False))
    print("\nblaVIM contig pairs:\n" + pairs.to_string(index=False))


if __name__ == "__main__":
    main()
