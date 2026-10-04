#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
pLIN v4.1 on short-read (Illumina-like) assemblies (exploratory, not pre-registered).

A  Fragmentation: 200 plasmids of the confirmatory test set (seed 5151). Reads are
   simulated with wgsim (2 x 150 bp, 60x, insert 400 bp, 0.2% sequencing error,
   no mutations; seed 7), assembled with SPAdes (--isolate), and the contigs of
   >= 500 bp are typed as one plasmid. Codes are compared level by level with the
   plasmid's published (complete-sequence) code.
B  End to end: the 8 Swiss VIM-1 genomes (long-read assemblies, chromosome and
   plasmids) -> simulated reads -> SPAdes -> MOB-recon plasmid bins -> pLIN v4.1
   (one session per run, all bins together).
C  As B with the real Illumina reads of the same 8 isolates (ENA PRJEB98563;
   Seth-Smith et al. 2026, Antimicrob Agents Chemother 70:e01827-25) instead of
   simulated reads. Each long-read plasmid is matched to
   the bin containing most of its k-mers, and their codes are compared, overall and
   separately for plasmids that MOB-recon reconstructed (>= 90% of their k-mers in
   one bin). L6 IDs new to the release are provisional and session-specific.

Usage:
  python v41_shortread.py [--index-dir DIR] [--part A|B|C|both|all]
Output: output/backbone_v41/shortread/{fragmentation_codes.tsv, fragmentation_metrics.tsv,
        vim_bins.tsv, vim_metrics.tsv, vim_real_bins.tsv, vim_real_metrics.tsv}
"""

import argparse
import glob
import os
import random
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from Bio import SeqIO

from plin_kmers import adaptive_sketch, containment
from plin_v41_typer import PlinV41Release

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "output", "backbone_v41", "shortread")
WGSIM = os.path.expanduser("~/miniconda3/envs/busco_env/bin/wgsim")
SPADES = os.path.expanduser("~/miniconda3/envs/mob_suite_env/bin/spades.py")
MOB_RECON = os.path.expanduser("~/miniconda3/envs/mob_suite_env/bin/mob_recon")
COVER, READ, INSERT, ERR = 60, 150, 400, 0.002


def simulate_and_assemble(name, fasta, workdir, threads=2):
    """wgsim reads -> SPAdes contigs (cached). Returns the contigs path or None."""
    d = os.path.join(workdir, name)
    contigs = os.path.join(d, "spades", "contigs.fasta")
    if os.path.exists(contigs):
        return contigs
    os.makedirs(d, exist_ok=True)
    length = sum(len(r.seq) for r in SeqIO.parse(fasta, "fasta"))
    pairs = max(1000, int(length * COVER / (2 * READ)))
    r1, r2 = os.path.join(d, "r1.fq"), os.path.join(d, "r2.fq")
    subprocess.run([WGSIM, "-S", "7", "-N", str(pairs), "-1", str(READ), "-2", str(READ), "-d", str(INSERT),
                    "-e", str(ERR), "-r", "0", "-R", "0", "-X", "0", fasta, r1, r2],
                   check=True, capture_output=True)
    r = subprocess.run([SPADES, "--isolate", "-1", r1, "-2", r2, "-o", os.path.join(d, "spades"), "-t", str(threads)],
                       capture_output=True)
    for f in (r1, r2):
        os.remove(f)
    return contigs if r.returncode == 0 and os.path.exists(contigs) else None


REAL_READS = os.path.join(BASE_DIR, "outbreak_validation", "swiss_vim1", "reads")


def assemble_real(name, workdir, threads=8):
    """Real Illumina reads (R1/R2) -> SPAdes contigs (cached). Returns the contigs path or None."""
    d = os.path.join(workdir, name)
    contigs = os.path.join(d, "spades", "contigs.fasta")
    if os.path.exists(contigs):
        return contigs
    os.makedirs(d, exist_ok=True)
    r1, r2 = (os.path.join(REAL_READS, f"{name}_R{i}.fastq.gz") for i in (1, 2))
    r = subprocess.run([SPADES, "--isolate", "-1", r1, "-2", r2, "-o", os.path.join(d, "spades"), "-t", str(threads)],
                       capture_output=True)
    return contigs if r.returncode == 0 and os.path.exists(contigs) else None


def shared(a, b):
    n = 0
    for x, y in zip(str(a).split("."), str(b).split(".")):
        if x != y:
            break
        n += 1
    return n


def part_a(rel, codes):
    test = pd.read_csv(os.path.join(BASE_DIR, "output", "backbone_v41", "confirm", "test_plasmids.tsv"), sep="\t")
    pick = sorted(random.Random(5151).sample(list(test.plasmid_id), 200))
    fasta = dict(zip(test.plasmid_id, test.fasta))
    work = os.path.join(OUT, "work_fragmentation")
    with ThreadPoolExecutor(6) as ex:
        contigs = dict(zip(pick, ex.map(lambda p: simulate_and_assemble(p, fasta[p], work), pick)))
    recs, info = [], {}
    for p in pick:
        cs = [str(r.seq) for r in SeqIO.parse(contigs[p], "fasta") if len(r.seq) >= 500] if contigs[p] else []
        length = sum(len(r.seq) for r in SeqIO.parse(fasta[p], "fasta"))
        info[p] = {"complete_length": length, "contigs": len(cs), "assembled_bp": sum(map(len, cs))}
        if cs:
            recs.append((p, ("N" * 50).join(cs)))
    typed = rel.type_sequences(recs, threads=8).set_index("plasmid_id")
    rows = []
    for p in pick:
        c = typed.pLIN_v41.get(p)
        rows.append({"plasmid_id": p, **info[p], "complete_code": codes[p], "shortread_code": c,
                     "levels_shared": shared(c, codes[p]) if c else 0})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "fragmentation_codes.tsv"), sep="\t", index=False)
    df["frag"] = pd.cut(df.contigs, [-1, 0, 1, 5, 20, 10 ** 6], labels=["not assembled", "1 contig", "2-5", "6-20", ">20"])
    met = pd.DataFrame([{"group": "all", "plasmids": len(df), **{f"L{k}_kept_pct": round(100 * (df.levels_shared >= k).mean(), 1) for k in range(1, 7)}}] +
                       [{"group": str(g), "plasmids": len(s), **{f"L{k}_kept_pct": round(100 * (s.levels_shared >= k).mean(), 1) for k in range(1, 7)}}
                        for g, s in df.groupby("frag", observed=True)])
    met.to_csv(os.path.join(OUT, "fragmentation_metrics.tsv"), sep="\t", index=False)
    print(met.to_string(index=False))


def part_b(rel, real=False):
    work = os.path.join(OUT, "work_vim_real" if real else "work_vim")
    tag = "vim_real" if real else "vim"
    asm = sorted(glob.glob(os.path.join(BASE_DIR, "outbreak_validation", "swiss_vim1", "assemblies", "*", "assembly.fasta")))
    lr = pd.read_csv(os.path.join(BASE_DIR, "output", "backbone_v41", "case_studies", "swiss_vim1_codes.tsv"), sep="\t")
    lr_seq = {}
    for f in glob.glob(os.path.join(BASE_DIR, "sample_data", "swiss_vim1_outbreak", "*_plasmids.fasta")):
        lr_seq.update({r.id: str(r.seq).upper() for r in SeqIO.parse(f, "fasta")})
    bins = []
    for a in asm:
        iso = os.path.basename(os.path.dirname(a))
        contigs = assemble_real(iso, work) if real else simulate_and_assemble(iso, a, work, threads=8)
        if contigs is None:
            raise RuntimeError(f"assembly failed for {iso}")
        mob = os.path.join(work, iso, "mob_recon")
        if not os.path.exists(os.path.join(mob, "contig_report.txt")):
            subprocess.run([MOB_RECON, "-i", contigs, "-o", mob, "-n", "8", "--force"], capture_output=True)
        for b in sorted(glob.glob(os.path.join(mob, "plasmid_*.fasta"))):
            seq = ("N" * 50).join(str(r.seq) for r in SeqIO.parse(b, "fasta"))
            bins.append((f"{iso}|{os.path.basename(b)[:-6]}", seq.upper()))
    typed = rel.type_sequences(bins, threads=8).set_index("plasmid_id")
    sk_bin = {b: adaptive_sketch(s) for b, s in bins}
    rows = []
    for _, r in lr.iterrows():
        iso = r.plasmid_id.split("_")[0]
        a_sk, a_sc = adaptive_sketch(lr_seq[r.plasmid_id])
        cand = [(containment(a_sk, *sk_bin[b], a_sc), b) for b, _ in bins if b.startswith(iso + "|")]
        best = max(cand) if cand else (0.0, None)
        c = typed.pLIN_v41.get(best[1]) if best[1] else None
        rows.append({"contig_id": r.plasmid_id, "length_bp": r.length_bp, "key_genes": r.key_genes,
                     "longread_code": r.pLIN_v41, "best_bin": best[1], "share_of_kmers_in_bin": round(best[0], 3),
                     "shortread_code": c, "levels_shared": shared(c, r.pLIN_v41) if c else 0})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, f"{tag}_bins.tsv"), sep="\t", index=False)
    vim = df[df.key_genes.fillna("").str.contains("blaVIM")]
    good = df[df.share_of_kmers_in_bin >= 0.9]                  # plasmid reconstructed by MOB-recon
    met = pd.DataFrame([{"set": name, "plasmids": len(s),
                         **{f"L{k}_kept_pct": round(100 * (s.levels_shared >= k).mean(), 1) for k in range(1, 7)}}
                        for name, s in (("all long-read plasmids", df), ("blaVIM-1 plasmids", vim),
                                        (">= 90% of k-mers in one MOB-recon bin", good),
                                        ("< 90% (binning split or merged the plasmid)", df[df.share_of_kmers_in_bin < 0.9]))])
    met.to_csv(os.path.join(OUT, f"{tag}_metrics.tsv"), sep="\t", index=False)
    print(df[["contig_id", "length_bp", "share_of_kmers_in_bin", "longread_code", "shortread_code", "levels_shared"]].to_string(index=False))
    print(met.to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index-dir", default=None)
    ap.add_argument("--part", default="both")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rel = PlinV41Release(os.path.join(BASE_DIR, "output", "backbone_v41", "release"), index_dir=args.index_dir)
    codes = pd.read_csv(os.path.join(BASE_DIR, "output", "backbone_v41", "codes_v41.tsv"), sep="\t") \
        .set_index("plasmid_id").pLIN_v41
    if args.part in ("A", "both", "all"):
        part_a(rel, codes)
    if args.part in ("B", "both", "all"):
        part_b(rel)
    if args.part in ("C", "all"):
        part_b(rel, real=True)


if __name__ == "__main__":
    main()
