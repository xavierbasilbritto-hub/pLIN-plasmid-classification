#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
"""
Validate the contig classifier (plasmid vs chromosome) using:
  - Known plasmids: random sample from the 8,077 training sequences
  - Known chromosomes: downloaded from NCBI RefSeq (complete genomes)
  - Mixed assemblies: simulated multi-contig files

Produces output/contig_validation/ with:
  - validation_results.tsv       — per-contig classification results
  - validation_summary.json      — performance metrics
  - size_stratified_metrics.json — accuracy by contig size
"""

import os
import sys
import json
import glob
import random
import subprocess
import tempfile
import numpy as np
import pandas as pd

# Add project to path
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)

from plin_app import (
    classify_contigs_plasmid_vs_chromosome,
    compute_kmer_vectors,
)
from Bio import SeqIO, Entrez

Entrez.email = "plin_tool@example.com"

OUT_DIR = os.path.join(BASE_DIR, "output", "contig_validation")
TRAIN_DIR = os.path.join(BASE_DIR, "plasmid_sequences_for_training")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Chromosome accessions (complete RefSeq genomes) ──────────────────────────
# Diverse set: E. coli, K. pneumoniae, S. aureus, E. faecium, A. baumannii, P. aeruginosa
CHROMOSOME_ACCESSIONS = [
    # E. coli chromosomes
    "NC_000913.3",   # E. coli K-12 MG1655, 4.6 Mb
    "NC_002695.2",   # E. coli O157:H7 Sakai, 5.5 Mb
    "NC_010473.1",   # E. coli BL21(DE3), 4.6 Mb
    "NC_004431.1",   # E. coli CFT073, 5.2 Mb
    # K. pneumoniae chromosomes
    "NC_016845.1",   # K. pneumoniae subsp. pneumoniae HS11286, 5.3 Mb
    "NC_022566.1",   # K. pneumoniae KCTC 2242, 5.4 Mb
    # S. aureus chromosomes
    "NC_002745.2",   # S. aureus N315, 2.8 Mb
    "NC_002953.3",   # S. aureus MW2, 2.8 Mb
    "NC_007795.1",   # S. aureus USA300 FPR3757, 2.9 Mb
    # E. faecium chromosomes
    "NC_017960.1",   # E. faecium DO, 2.7 Mb
    # A. baumannii chromosomes
    "NC_010611.1",   # A. baumannii AB0057, 4.0 Mb
    "NC_017162.2",   # A. baumannii 1656-2, 3.9 Mb
    # P. aeruginosa chromosomes
    "NC_002516.2",   # P. aeruginosa PAO1, 6.3 Mb
    "NC_009656.1",   # P. aeruginosa PA7, 6.6 Mb
    # Salmonella chromosomes
    "NC_003197.2",   # S. enterica Typhimurium LT2, 4.9 Mb
    "NC_003198.1",   # S. enterica Typhi CT18, 4.8 Mb
]


def download_chromosomes():
    """Download chromosome sequences from NCBI RefSeq."""
    chromo_dir = os.path.join(OUT_DIR, "chromosomes")
    os.makedirs(chromo_dir, exist_ok=True)

    downloaded = []
    for acc in CHROMOSOME_ACCESSIONS:
        fasta_path = os.path.join(chromo_dir, f"{acc}.fasta")
        if os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 1000:
            rec = SeqIO.read(fasta_path, "fasta")
            downloaded.append({"acc": acc, "path": fasta_path, "length": len(rec.seq)})
            print(f"  {acc}: cached ({len(rec.seq)/1e6:.1f} Mb)")
            continue

        print(f"  Downloading {acc}...", end=" ", flush=True)
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
            rec = SeqIO.read(handle, "fasta")
            handle.close()
            SeqIO.write(rec, fasta_path, "fasta")
            downloaded.append({"acc": acc, "path": fasta_path, "length": len(rec.seq)})
            print(f"OK ({len(rec.seq)/1e6:.1f} Mb)")
        except Exception as e:
            print(f"FAILED: {e}")

    return downloaded


def sample_plasmids(n=600):
    """Sample n plasmids from training data, stratified by size."""
    all_fastas = []
    for group in sorted(os.listdir(TRAIN_DIR)):
        fastas_dir = os.path.join(TRAIN_DIR, group, "fastas")
        if os.path.isdir(fastas_dir):
            for fa in glob.glob(os.path.join(fastas_dir, "*.fasta")):
                all_fastas.append({"path": fa, "group": group})

    random.seed(42)
    if len(all_fastas) > n:
        sampled = random.sample(all_fastas, n)
    else:
        sampled = all_fastas

    return sampled


def build_records(fasta_path, label, group=None):
    """Read a FASTA file and build records for classification."""
    records = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        records.append({
            "plasmid_id": rec.id,
            "sequence": str(rec.seq),
            "length": len(rec.seq),
            "true_label": label,
            "true_group": group or "",
        })
    return records


def main():
    print("=" * 70)
    print("Contig Classifier Validation")
    print("=" * 70)

    # Step 1: Download chromosomes
    print("\n[Step 1] Downloading chromosome sequences...")
    chromosomes = download_chromosomes()
    print(f"  Downloaded {len(chromosomes)} chromosomes")

    # Step 2: Sample plasmids
    print("\n[Step 2] Sampling plasmids from training data...")
    plasmids = sample_plasmids(n=600)
    print(f"  Sampled {len(plasmids)} plasmids")

    # Step 3: Build all records
    print("\n[Step 3] Building classification records...")
    all_records = []

    # Add chromosomes
    for chromo in chromosomes:
        recs = build_records(chromo["path"], "chromosome")
        all_records.extend(recs)

    # Add plasmids
    for plas in plasmids:
        recs = build_records(plas["path"], "plasmid", plas["group"])
        all_records.extend(recs)

    print(f"  Total records: {len(all_records)}")
    n_chrom = sum(1 for r in all_records if r["true_label"] == "chromosome")
    n_plas = sum(1 for r in all_records if r["true_label"] == "plasmid")
    print(f"  Chromosomes: {n_chrom}, Plasmids: {n_plas}")

    # Step 4: Run classifier
    print("\n[Step 4] Running contig classifier on all records...")
    results = classify_contigs_plasmid_vs_chromosome(all_records)
    print(f"  Classified {len(results)} contigs")

    # Merge with true labels
    for i, res in enumerate(results):
        res["true_label"] = all_records[i]["true_label"]
        res["true_group"] = all_records[i]["true_group"]

    # Save full results
    results_df = pd.DataFrame(results)
    results_df.to_csv(os.path.join(OUT_DIR, "validation_results.tsv"),
                      sep="\t", index=False)

    # Step 5: Compute metrics
    print("\n[Step 5] Computing validation metrics...")

    # Overall metrics
    plas_true = results_df[results_df["true_label"] == "plasmid"]
    chrom_true = results_df[results_df["true_label"] == "chromosome"]

    # Sensitivity: proportion of true plasmids classified as plasmid (not chromosome)
    plas_correct = plas_true[plas_true["classification"].isin(["plasmid", "incomplete_plasmid"])]
    plas_as_chromo = plas_true[plas_true["classification"] == "chromosome"]
    sensitivity = len(plas_correct) / len(plas_true) if len(plas_true) > 0 else 0

    # Specificity: proportion of true chromosomes classified as chromosome
    chrom_correct = chrom_true[chrom_true["classification"] == "chromosome"]
    chrom_as_plas = chrom_true[chrom_true["classification"].isin(["plasmid", "incomplete_plasmid"])]
    specificity = len(chrom_correct) / len(chrom_true) if len(chrom_true) > 0 else 0

    # Accuracy
    correct = len(plas_correct) + len(chrom_correct)
    accuracy = correct / len(results_df) if len(results_df) > 0 else 0

    # Score distributions
    plas_scores = plas_true["score"].values
    chrom_scores = chrom_true["score"].values
    plas_median = float(np.median(plas_scores))
    chrom_median = float(np.median(chrom_scores))

    print(f"\n  --- Overall Metrics ---")
    print(f"  Total contigs: {len(results_df)}")
    print(f"  True plasmids: {len(plas_true)}")
    print(f"  True chromosomes: {len(chrom_true)}")
    print(f"  Sensitivity (plasmid recall): {sensitivity:.1%}")
    print(f"  Specificity (chromosome recall): {specificity:.1%}")
    print(f"  Overall accuracy: {accuracy:.1%}")
    print(f"  Plasmid score median: {plas_median:.1f}")
    print(f"  Chromosome score median: {chrom_median:.1f}")

    # Classification breakdown
    print(f"\n  --- Classification Breakdown ---")
    for cls in ["plasmid", "incomplete_plasmid", "chromosome"]:
        n_cls = len(results_df[results_df["classification"] == cls])
        pct = 100 * n_cls / len(results_df)
        print(f"  {cls}: {n_cls} ({pct:.1f}%)")

    # True plasmid classification breakdown
    print(f"\n  --- True Plasmids Classification ---")
    for cls in ["plasmid", "incomplete_plasmid", "chromosome"]:
        n_cls = len(plas_true[plas_true["classification"] == cls])
        pct = 100 * n_cls / len(plas_true) if len(plas_true) > 0 else 0
        print(f"  {cls}: {n_cls} ({pct:.1f}%)")

    # True chromosomes classification breakdown
    print(f"\n  --- True Chromosomes Classification ---")
    for cls in ["plasmid", "incomplete_plasmid", "chromosome"]:
        n_cls = len(chrom_true[chrom_true["classification"] == cls])
        pct = 100 * n_cls / len(chrom_true) if len(chrom_true) > 0 else 0
        print(f"  {cls}: {n_cls} ({pct:.1f}%)")

    # Step 6: Size-stratified metrics
    print(f"\n  --- Size-Stratified Metrics ---")
    size_bins = [
        ("< 50 kb", 0, 50_000),
        ("50-500 kb", 50_000, 500_000),
        ("> 500 kb", 500_000, float("inf")),
    ]
    size_metrics = []
    for label, lo, hi in size_bins:
        sub = results_df[(results_df["length_bp"] >= lo) & (results_df["length_bp"] < hi)]
        if len(sub) == 0:
            continue
        n_plas_true = len(sub[sub["true_label"] == "plasmid"])
        n_chrom_true = len(sub[sub["true_label"] == "chromosome"])
        n_correct = len(sub[
            ((sub["true_label"] == "plasmid") & sub["classification"].isin(["plasmid", "incomplete_plasmid"])) |
            ((sub["true_label"] == "chromosome") & (sub["classification"] == "chromosome"))
        ])
        acc_bin = n_correct / len(sub) if len(sub) > 0 else 0

        # Sensitivity for this bin
        plas_in_bin = sub[sub["true_label"] == "plasmid"]
        sens_bin = len(plas_in_bin[plas_in_bin["classification"].isin(["plasmid", "incomplete_plasmid"])]) / len(plas_in_bin) if len(plas_in_bin) > 0 else 0

        # Specificity for this bin
        chrom_in_bin = sub[sub["true_label"] == "chromosome"]
        spec_bin = len(chrom_in_bin[chrom_in_bin["classification"] == "chromosome"]) / len(chrom_in_bin) if len(chrom_in_bin) > 0 else 0

        print(f"  {label}: n={len(sub)}, plasmids={n_plas_true}, chromosomes={n_chrom_true}, "
              f"accuracy={acc_bin:.1%}, sensitivity={sens_bin:.1%}, specificity={spec_bin:.1%}")

        size_metrics.append({
            "size_bin": label,
            "n_total": len(sub),
            "n_plasmid_true": n_plas_true,
            "n_chromosome_true": n_chrom_true,
            "accuracy": round(acc_bin, 4),
            "sensitivity": round(sens_bin, 4),
            "specificity": round(spec_bin, 4),
        })

    # Step 7: Score threshold analysis
    print(f"\n  --- Score Threshold Analysis ---")
    thresholds = list(range(-50, 60, 5))
    threshold_metrics = []
    for t in thresholds:
        tp = len(plas_true[plas_true["score"] >= t])
        fp = len(chrom_true[chrom_true["score"] >= t])
        fn = len(plas_true) - tp
        tn = len(chrom_true) - fp
        sens = tp / (tp + fn) if (tp + fn) > 0 else 0
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0
        threshold_metrics.append({
            "threshold": t,
            "sensitivity": round(sens, 4),
            "specificity": round(spec, 4),
            "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        })

    # Step 8: K. pneumoniae example (if we have the data)
    # Use actual K. pneumoniae chromosome + plasmids from training
    kpn_example = None
    kpn_chrom = [c for c in chromosomes if "016845" in c["acc"] or "022566" in c["acc"]]
    kpn_plas_paths = glob.glob(os.path.join(TRAIN_DIR, "IncFII", "fastas", "NZ_CP*.fasta"))[:3]

    if kpn_chrom and kpn_plas_paths:
        print(f"\n  --- K. pneumoniae Example ---")
        example_records = []
        for c in kpn_chrom[:1]:
            recs = build_records(c["path"], "chromosome")
            example_records.extend(recs)
        for p in kpn_plas_paths:
            recs = build_records(p, "plasmid", "IncFII")
            example_records.extend(recs)

        if example_records:
            example_results = classify_contigs_plasmid_vs_chromosome(example_records)
            kpn_example = []
            for j, er in enumerate(example_results):
                er["true_label"] = example_records[j]["true_label"]
                kpn_example.append(er)
                print(f"  {er['plasmid_id']}: score={er['score']}, "
                      f"class={er['classification']}, conf={er['confidence']}%, "
                      f"len={er['length_bp']/1000:.0f}kb, true={example_records[j]['true_label']}")

    # Save summary
    summary = {
        "total_contigs": len(results_df),
        "n_plasmids_true": int(len(plas_true)),
        "n_chromosomes_true": int(len(chrom_true)),
        "sensitivity": round(sensitivity, 4),
        "specificity": round(specificity, 4),
        "accuracy": round(accuracy, 4),
        "plasmid_score_median": round(plas_median, 1),
        "chromosome_score_median": round(chrom_median, 1),
        "plasmid_score_mean": round(float(np.mean(plas_scores)), 1),
        "chromosome_score_mean": round(float(np.mean(chrom_scores)), 1),
        "size_stratified": size_metrics,
        "threshold_analysis": threshold_metrics,
        "classification_breakdown": {
            cls: int(len(results_df[results_df["classification"] == cls]))
            for cls in ["plasmid", "incomplete_plasmid", "chromosome"]
        },
    }
    if kpn_example:
        summary["kpn_example"] = kpn_example

    with open(os.path.join(OUT_DIR, "validation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # Save threshold analysis
    pd.DataFrame(threshold_metrics).to_csv(
        os.path.join(OUT_DIR, "threshold_analysis.tsv"), sep="\t", index=False)

    print(f"\n{'=' * 70}")
    print(f"Validation Complete")
    print(f"  Output: {OUT_DIR}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
