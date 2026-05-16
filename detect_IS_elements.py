#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
"""
Detect IS elements in pLIN training plasmids using blastn against
curated IS element reference sequences from ISfinder/NCBI.

Outputs:
  output/mge_detection/is_element_hits.tsv       — all IS hits per plasmid
  output/mge_detection/is_family_counts.tsv      — IS family × organism category
  output/mge_detection/composite_transposons.tsv — IS pairs flanking AMR genes
  output/mge_detection/is_summary.json           — summary statistics
"""

import os
import sys
import json
import subprocess
import tempfile
import glob
import re
from collections import defaultdict, Counter

import pandas as pd
import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TRAIN_DIR = os.path.join(BASE_DIR, "plasmid_sequences_for_training")
AMR_FILE = os.path.join(BASE_DIR, "output", "amrfinder", "amrfinder_all_plasmids.tsv")
OUT_DIR = os.path.join(BASE_DIR, "output", "mge_detection")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Curated IS element reference sequences ──────────────────────────────────
# Representative transposase sequences for the major IS families found in
# Enterobacteriaceae, Staphylococcus, and Enterococcus plasmids.
# Accessions from NCBI nucleotide database (verified against ISfinder).

IS_REFERENCES = {
    # Gram-negative IS families
    "IS26":    {"acc": "X00011.1", "len_range": (820, 820),   "desc": "IS26 transposase"},
    "ISEcp1":  {"acc": "AJ242809.1", "len_range": (1656, 1656), "desc": "ISEcp1 transposase"},
    "IS1":     {"acc": "J01730.1", "len_range": (768, 768),   "desc": "IS1 transposase"},
    "IS903":   {"acc": "M17148.1", "len_range": (1057, 1057), "desc": "IS903 transposase"},
    "IS6100":  {"acc": "M95400.1", "len_range": (880, 880),   "desc": "IS6100 transposase"},
    "ISKpn26": {"acc": "KF914891.1", "len_range": (1200, 1500), "desc": "ISKpn26 family"},
    "IS5":     {"acc": "X02311.1", "len_range": (1195, 1195), "desc": "IS5 transposase"},
    "IS3":     {"acc": "X02180.1", "len_range": (1258, 1258), "desc": "IS3 transposase"},
    "IS4321":  {"acc": "AJ245418.1", "len_range": (1700, 1700), "desc": "IS4321 transposase"},
    "IS15":    {"acc": "X01840.1", "len_range": (1400, 1400), "desc": "IS15DI transposase"},
    "IS10":    {"acc": "J01830.1", "len_range": (1329, 1329), "desc": "IS10R transposase"},
    "IS2":     {"acc": "J01733.1", "len_range": (1327, 1327), "desc": "IS2 transposase"},
    "IS4":     {"acc": "V00029.1", "len_range": (1426, 1426), "desc": "IS4 transposase"},
    "IS30":    {"acc": "X00792.1", "len_range": (1221, 1221), "desc": "IS30 transposase"},
    "IS66":    {"acc": "X53365.1", "len_range": (2548, 2548), "desc": "IS66 transposase"},
    "IS110":   {"acc": "M21395.1", "len_range": (1449, 1449), "desc": "IS110 transposase"},
    "ISPa":    {"acc": "AF261825.1", "len_range": (1000, 2000), "desc": "ISPa family (P. aeruginosa)"},
    "ISAba":   {"acc": "AY758396.1", "len_range": (1000, 2000), "desc": "ISAba family (A. baumannii)"},
    # Gram-positive IS families
    "IS256":   {"acc": "M18086.1", "len_range": (1324, 1324), "desc": "IS256 transposase"},
    "IS257":   {"acc": "U40412.1", "len_range": (789, 789),   "desc": "IS257/IS431 transposase"},
    "IS16":    {"acc": "AF053365.1", "len_range": (1300, 1300), "desc": "IS16 transposase (Enterococcus)"},
    "ISEnfa":  {"acc": "AF162694.1", "len_range": (1200, 1200), "desc": "ISEnfa transposase (E. faecalis)"},
    "IS1216":  {"acc": "L40841.1", "len_range": (813, 813),   "desc": "IS1216V transposase"},
    "IS1251":  {"acc": "X83579.1", "len_range": (1300, 1300), "desc": "IS1251 (S. aureus)"},
    "Tn916":   {"acc": "U09422.1", "len_range": (18000, 18000), "desc": "Tn916 conjugative transposon"},
}

# Organism categories for each training group
GRAM_CATEGORY = {
    "ColE": "Gram-negative", "ColRNAI": "Gram-negative", "IncA": "Gram-negative",
    "IncAC2": "Gram-negative", "IncC": "Gram-negative", "IncF": "Gram-negative",
    "IncFIB": "Gram-negative", "IncFIBK": "Gram-negative", "IncFIC": "Gram-negative",
    "IncFII": "Gram-negative", "IncHI1": "Gram-negative", "IncHI2": "Gram-negative",
    "IncI": "Gram-negative", "IncI1": "Gram-negative", "IncI2": "Gram-negative",
    "IncN": "Gram-negative", "IncR": "Gram-negative", "IncX1": "Gram-negative",
    "IncX3": "Gram-negative", "IncX4": "Gram-negative",
    "repSA_large": "Gram-positive", "repSA_small": "Gram-positive",
    "repEF_conj": "Gram-positive", "repEF_res": "Gram-positive",
    "repAci1": "Acinetobacter", "repAci_large": "Acinetobacter",
    "repPae_large": "Pseudomonas", "repPae_small": "Pseudomonas",
}


def download_is_references(out_fasta):
    """Download IS element reference sequences from NCBI using efetch."""
    from Bio import Entrez, SeqIO
    Entrez.email = "plin_tool@example.com"

    records = []
    for is_name, info in IS_REFERENCES.items():
        acc = info["acc"]
        print(f"  Fetching {is_name} ({acc})...", end=" ", flush=True)
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
            rec = SeqIO.read(handle, "fasta")
            handle.close()
            # Rename header to IS family name
            rec.id = is_name
            rec.description = info["desc"]
            records.append(rec)
            print(f"OK ({len(rec.seq)} bp)")
        except Exception as e:
            print(f"FAILED: {e}")

    with open(out_fasta, "w") as f:
        SeqIO.write(records, f, "fasta")
    print(f"\nSaved {len(records)} IS reference sequences to {out_fasta}")
    return len(records)


def build_blast_db(fasta_path):
    """Build a BLAST nucleotide database from the IS reference FASTA."""
    cmd = ["makeblastdb", "-in", fasta_path, "-dbtype", "nucl",
           "-parse_seqids", "-out", fasta_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"makeblastdb error: {result.stderr}")
        sys.exit(1)
    print("BLAST database built successfully.")


def scan_plasmids_for_IS(is_db_path, training_dir, out_tsv):
    """Run blastn on all training plasmids against the IS reference database.

    Uses megablast with relaxed parameters suitable for IS element detection:
    - evalue 1e-10: stringent enough to avoid false positives
    - perc_identity 80: IS elements can have 80-100% identity across species
    - qcov_hsp_perc 50: at least half the query IS element must align
    """
    all_hits = []
    groups = sorted(os.listdir(training_dir))
    total_plasmids = 0
    plasmids_with_IS = 0

    for group in groups:
        fastas_dir = os.path.join(training_dir, group, "fastas")
        if not os.path.isdir(fastas_dir):
            continue

        fasta_files = sorted(glob.glob(os.path.join(fastas_dir, "*.fasta")))
        if not fasta_files:
            continue

        print(f"\n  Scanning {group} ({len(fasta_files)} plasmids)...", flush=True)

        # Concatenate all FASTAs for this group into a temp file for batch BLAST
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
            tmp_path = tmp.name
            for fa in fasta_files:
                with open(fa) as f:
                    for line in f:
                        tmp.write(line)

        # Run blastn
        cmd = [
            "blastn",
            "-query", tmp_path,
            "-db", is_db_path,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
            "-evalue", "1e-10",
            "-perc_identity", "80",
            "-max_target_seqs", "50",
            "-num_threads", "4",
            "-task", "blastn",  # more sensitive than megablast for IS detection
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        os.unlink(tmp_path)

        if result.returncode != 0:
            print(f"    BLAST error for {group}: {result.stderr[:200]}")
            continue

        # Parse BLAST output
        group_hits = 0
        group_plasmids_with_hits = set()
        for line in result.stdout.strip().split("\n"):
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 14:
                continue

            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen = parts

            # Filter: alignment must cover ≥50% of the IS reference
            align_len = int(length)
            ref_len = int(slen)
            coverage = align_len / ref_len if ref_len > 0 else 0

            if coverage >= 0.50 and float(pident) >= 80.0:
                # Map accession ID to source file name
                source_file = qseqid
                is_family = sseqid

                all_hits.append({
                    "source_file": source_file,
                    "inc_group": group,
                    "organism_category": GRAM_CATEGORY.get(group, "Unknown"),
                    "is_family": is_family,
                    "pident": float(pident),
                    "align_length": align_len,
                    "ref_length": ref_len,
                    "coverage": round(coverage, 3),
                    "q_start": int(qstart),
                    "q_end": int(qend),
                    "s_start": int(sstart),
                    "s_end": int(send),
                    "evalue": float(evalue),
                    "bitscore": float(bitscore),
                    "plasmid_length": int(qlen),
                })
                group_hits += 1
                group_plasmids_with_hits.add(source_file)

        total_plasmids += len(fasta_files)
        plasmids_with_IS += len(group_plasmids_with_hits)
        print(f"    {group}: {group_hits} IS hits in {len(group_plasmids_with_hits)}/{len(fasta_files)} plasmids")

    # Save all hits
    df = pd.DataFrame(all_hits)
    if len(df) > 0:
        df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\n  Total: {len(df)} IS hits across {plasmids_with_IS}/{total_plasmids} plasmids")
    return df


def detect_composite_transposons(is_hits_df, amr_df):
    """Identify composite transposons: paired IS elements flanking AMR genes.

    Logic: For each plasmid, find IS element pairs of the same family where
    an AMR gene falls between them (within coordinates).
    """
    if is_hits_df is None or len(is_hits_df) == 0:
        return pd.DataFrame()

    # Load AMR gene positions per plasmid
    amr_positions = {}
    if amr_df is not None and len(amr_df) > 0:
        for _, row in amr_df.iterrows():
            source = str(row.get("source_file", ""))
            # Extract accession from source_file (strip prefix)
            acc = source.split("_", 1)[-1] if "_" in source else source
            start = int(row.get("Start", 0))
            end = int(row.get("Stop", 0))
            symbol = str(row.get("Element symbol", ""))
            gene_type = str(row.get("Type", ""))
            if gene_type == "AMR" and start > 0:
                if acc not in amr_positions:
                    amr_positions[acc] = []
                amr_positions[acc].append({"start": start, "end": end, "gene": symbol})

    composite_transposons = []

    for source_file, group in is_hits_df.groupby("source_file"):
        # Get all IS hits for this plasmid, sorted by position
        hits = group.sort_values("q_start")
        is_families_present = hits["is_family"].unique()

        for is_fam in is_families_present:
            fam_hits = hits[hits["is_family"] == is_fam].reset_index(drop=True)
            if len(fam_hits) < 2:
                continue

            # Check each pair of same-family IS elements
            for i in range(len(fam_hits) - 1):
                is1_end = fam_hits.iloc[i]["q_end"]
                is2_start = fam_hits.iloc[i + 1]["q_start"]

                # Check if any AMR gene falls between these IS elements
                acc = source_file
                if acc in amr_positions:
                    cargo_genes = [g for g in amr_positions[acc]
                                   if g["start"] >= is1_end and g["end"] <= is2_start]
                    if cargo_genes:
                        composite_transposons.append({
                            "source_file": source_file,
                            "inc_group": fam_hits.iloc[i]["inc_group"],
                            "organism_category": fam_hits.iloc[i]["organism_category"],
                            "is_family": is_fam,
                            "is1_start": fam_hits.iloc[i]["q_start"],
                            "is1_end": int(is1_end),
                            "is2_start": int(is2_start),
                            "is2_end": fam_hits.iloc[i + 1]["q_end"],
                            "cargo_genes": ",".join(g["gene"] for g in cargo_genes),
                            "n_cargo_genes": len(cargo_genes),
                            "region_size": int(is2_start) - int(is1_end),
                        })

    ct_df = pd.DataFrame(composite_transposons)
    return ct_df


def main():
    print("=" * 70)
    print("IS Element Detection in pLIN Training Plasmids")
    print("=" * 70)

    # Step 1: Download IS reference sequences
    is_ref_fasta = os.path.join(OUT_DIR, "is_reference_sequences.fasta")
    if not os.path.exists(is_ref_fasta) or os.path.getsize(is_ref_fasta) < 100:
        print("\n[Step 1] Downloading IS element reference sequences from NCBI...")
        n_refs = download_is_references(is_ref_fasta)
        if n_refs == 0:
            print("ERROR: No IS references downloaded. Check internet connection.")
            sys.exit(1)
    else:
        print(f"\n[Step 1] Using existing IS references: {is_ref_fasta}")

    # Step 2: Build BLAST database
    print("\n[Step 2] Building BLAST database...")
    build_blast_db(is_ref_fasta)

    # Step 3: Scan all training plasmids
    print("\n[Step 3] Scanning 8,077 training plasmids for IS elements...")
    hits_tsv = os.path.join(OUT_DIR, "is_element_hits.tsv")
    is_df = scan_plasmids_for_IS(is_ref_fasta, TRAIN_DIR, hits_tsv)

    if len(is_df) == 0:
        print("WARNING: No IS elements detected. Check BLAST database and parameters.")
        return

    # Step 4: Summarize IS family counts
    print("\n[Step 4] Summarizing IS family counts...")

    # Per IS family × organism category
    family_counts = is_df.groupby(["is_family", "organism_category"]).agg(
        n_hits=("source_file", "count"),
        n_plasmids=("source_file", "nunique"),
    ).reset_index()
    family_counts.to_csv(os.path.join(OUT_DIR, "is_family_counts.tsv"),
                         sep="\t", index=False)

    # Top 15 IS families overall
    top15 = is_df.groupby("is_family").agg(
        total_hits=("source_file", "count"),
        n_plasmids=("source_file", "nunique"),
    ).sort_values("total_hits", ascending=False).head(15)
    print("\nTop 15 IS families:")
    print(top15.to_string())

    # Per Inc/rep group
    group_counts = is_df.groupby("inc_group").agg(
        total_IS_hits=("source_file", "count"),
        n_plasmids_with_IS=("source_file", "nunique"),
    ).sort_values("total_IS_hits", ascending=False)
    group_counts.to_csv(os.path.join(OUT_DIR, "is_by_inc_group.tsv"),
                        sep="\t", index=False)

    # Step 5: Detect composite transposons
    print("\n[Step 5] Detecting composite transposons (IS pairs flanking AMR genes)...")
    amr_df = None
    if os.path.exists(AMR_FILE):
        amr_df = pd.read_csv(AMR_FILE, sep="\t")
        print(f"  Loaded {len(amr_df)} AMR gene annotations")

    ct_df = detect_composite_transposons(is_df, amr_df)
    if len(ct_df) > 0:
        ct_df.to_csv(os.path.join(OUT_DIR, "composite_transposons.tsv"),
                     sep="\t", index=False)
        print(f"\n  Found {len(ct_df)} composite transposons")

        # Top IS-AMR associations
        if "cargo_genes" in ct_df.columns:
            is_amr_pairs = []
            for _, row in ct_df.iterrows():
                for gene in row["cargo_genes"].split(","):
                    is_amr_pairs.append(f"{row['is_family']}-{gene}")
            pair_counts = Counter(is_amr_pairs).most_common(10)
            print("\n  Top IS-AMR associations:")
            for pair, count in pair_counts:
                print(f"    {pair}: n={count}")

        # CT by group
        ct_by_group = ct_df.groupby("inc_group").size().sort_values(ascending=False)
        print("\n  Composite transposons by Inc/rep group:")
        print(ct_by_group.to_string())
    else:
        print("  No composite transposons detected.")

    # Step 6: Save summary JSON
    total_plasmids = 8077
    plasmids_with_IS = is_df["source_file"].nunique()

    summary = {
        "total_training_plasmids": total_plasmids,
        "plasmids_with_IS_elements": plasmids_with_IS,
        "total_IS_hits": len(is_df),
        "unique_IS_families_detected": int(is_df["is_family"].nunique()),
        "top_15_IS_families": top15.reset_index().to_dict(orient="records"),
        "composite_transposons_total": len(ct_df) if len(ct_df) > 0 else 0,
        "sensitivity_note": "Detection based on blastn against curated IS reference sequences (≥80% identity, ≥50% coverage)",
    }
    with open(os.path.join(OUT_DIR, "is_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\n{'=' * 70}")
    print(f"IS Detection Complete")
    print(f"  Plasmids scanned: {total_plasmids}")
    print(f"  Plasmids with IS: {plasmids_with_IS} ({100*plasmids_with_IS/total_plasmids:.1f}%)")
    print(f"  Total IS hits: {len(is_df)}")
    print(f"  IS families found: {is_df['is_family'].nunique()}")
    print(f"  Composite transposons: {len(ct_df) if len(ct_df) > 0 else 0}")
    print(f"  Output directory: {OUT_DIR}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
