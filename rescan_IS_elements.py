#!/usr/bin/env python3
"""
Re-scan training plasmids with curated IS reference database.
Uses stricter parameters to avoid false positives from oversized references.
"""

import os
import sys
import json
import subprocess
import tempfile
import glob
from collections import Counter

import pandas as pd
import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TRAIN_DIR = os.path.join(BASE_DIR, "plasmid_sequences_for_training")
AMR_FILE = os.path.join(BASE_DIR, "output", "amrfinder", "amrfinder_all_plasmids.tsv")
OUT_DIR = os.path.join(BASE_DIR, "output", "mge_detection")
IS_DB = os.path.join(OUT_DIR, "is_reference_combined.fasta")

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


def scan_group(group, fastas_dir, is_db):
    """Scan all plasmids in a group for IS elements."""
    fasta_files = sorted(glob.glob(os.path.join(fastas_dir, "*.fasta")))
    if not fasta_files:
        return [], 0

    # Concatenate all FASTAs for batch BLAST
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        tmp_path = tmp.name
        for fa in fasta_files:
            with open(fa) as f:
                for line in f:
                    tmp.write(line)

    cmd = [
        "blastn",
        "-query", tmp_path,
        "-db", is_db,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
        "-evalue", "1e-20",        # stricter e-value
        "-perc_identity", "85",     # stricter identity
        "-max_target_seqs", "50",
        "-num_threads", "4",
        "-task", "blastn",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    os.unlink(tmp_path)

    if result.returncode != 0:
        print(f"    BLAST error: {result.stderr[:200]}")
        return [], len(fasta_files)

    hits = []
    for line in result.stdout.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 14:
            continue

        qseqid, sseqid, pident, length, mm, go, qstart, qend, sstart, send, evalue, bitscore, qlen, slen = parts

        align_len = int(length)
        ref_len = int(slen)
        coverage = align_len / ref_len if ref_len > 0 else 0

        # Require ≥70% coverage of the IS reference AND ≥85% identity
        if coverage >= 0.70 and float(pident) >= 85.0:
            hits.append({
                "source_file": qseqid,
                "inc_group": group,
                "organism_category": GRAM_CATEGORY.get(group, "Unknown"),
                "is_family": sseqid,
                "pident": float(pident),
                "align_length": align_len,
                "ref_length": ref_len,
                "coverage": round(coverage, 3),
                "q_start": int(qstart),
                "q_end": int(qend),
                "evalue": float(evalue),
                "bitscore": float(bitscore),
                "plasmid_length": int(qlen),
            })

    return hits, len(fasta_files)


def detect_composite_transposons(is_df, amr_file):
    """Find IS pairs flanking AMR genes."""
    amr_df = pd.read_csv(amr_file, sep="\t")
    amr_only = amr_df[amr_df["Type"] == "AMR"].copy()

    # Build AMR position lookup: source_file -> list of (start, end, gene)
    amr_positions = {}
    for _, row in amr_only.iterrows():
        sf = str(row["source_file"])
        # The AMR file source_file format: "RefSeq_NZ_KX881941.1" or similar
        # The BLAST query ID format: "RefSeq_NZ_KX881941.1" (from FASTA headers)
        acc = sf
        start = int(row.get("Start", 0))
        end = int(row.get("Stop", 0))
        gene = str(row.get("Element symbol", ""))
        if start > 0:
            amr_positions.setdefault(acc, []).append({"start": start, "end": end, "gene": gene})

    # Also try with just the accession part (strip "RefSeq_" prefix)
    for key in list(amr_positions.keys()):
        if key.startswith("RefSeq_"):
            short = key[7:]  # strip "RefSeq_"
            if short not in amr_positions:
                amr_positions[short] = amr_positions[key]

    ct_results = []
    for source_file, group in is_df.groupby("source_file"):
        hits = group.sort_values("q_start")
        for is_fam in hits["is_family"].unique():
            fam_hits = hits[hits["is_family"] == is_fam].reset_index(drop=True)
            if len(fam_hits) < 2:
                continue

            for i in range(len(fam_hits) - 1):
                is1_end = max(fam_hits.iloc[i]["q_start"], fam_hits.iloc[i]["q_end"])
                is2_start = min(fam_hits.iloc[i+1]["q_start"], fam_hits.iloc[i+1]["q_end"])

                # Region between IS elements shouldn't be too large (< 50kb)
                region = abs(is2_start - is1_end)
                if region > 50000 or region < 100:
                    continue

                # Check for AMR genes in between
                for key_variant in [source_file, f"RefSeq_{source_file}", source_file.replace("RefSeq_", "")]:
                    if key_variant in amr_positions:
                        cargo = [g for g in amr_positions[key_variant]
                                 if g["start"] >= is1_end and g["end"] <= is2_start]
                        if cargo:
                            ct_results.append({
                                "source_file": source_file,
                                "inc_group": fam_hits.iloc[i]["inc_group"],
                                "organism_category": fam_hits.iloc[i]["organism_category"],
                                "is_family": is_fam,
                                "is1_pos": f"{fam_hits.iloc[i]['q_start']}-{fam_hits.iloc[i]['q_end']}",
                                "is2_pos": f"{fam_hits.iloc[i+1]['q_start']}-{fam_hits.iloc[i+1]['q_end']}",
                                "region_size": region,
                                "cargo_genes": ",".join(g["gene"] for g in cargo),
                                "n_cargo": len(cargo),
                            })
                            break

    return pd.DataFrame(ct_results)


def main():
    print("=" * 70)
    print("IS Element Re-scan with Curated Reference Database")
    print("=" * 70)

    all_hits = []
    total_plasmids = 0
    total_with_IS = 0

    groups = sorted(os.listdir(TRAIN_DIR))
    for group in groups:
        fastas_dir = os.path.join(TRAIN_DIR, group, "fastas")
        if not os.path.isdir(fastas_dir):
            continue

        print(f"\n  {group}...", end=" ", flush=True)
        hits, n_plasmids = scan_group(group, fastas_dir, IS_DB)
        total_plasmids += n_plasmids

        n_with_IS = len(set(h["source_file"] for h in hits))
        total_with_IS += n_with_IS
        print(f"{len(hits)} hits in {n_with_IS}/{n_plasmids} plasmids")
        all_hits.extend(hits)

    is_df = pd.DataFrame(all_hits)
    is_df.to_csv(os.path.join(OUT_DIR, "is_element_hits_curated.tsv"), sep="\t", index=False)

    print(f"\n{'=' * 70}")
    print(f"RESULTS SUMMARY")
    print(f"{'=' * 70}")
    print(f"Total plasmids scanned: {total_plasmids}")
    print(f"Plasmids with IS elements: {total_with_IS} ({100*total_with_IS/total_plasmids:.1f}%)")
    print(f"Total IS hits: {len(is_df)}")
    print(f"Unique IS families: {is_df['is_family'].nunique()}")

    # Top IS families
    print(f"\n--- Top IS families ---")
    fam_counts = is_df.groupby("is_family").agg(
        total_hits=("source_file", "count"),
        unique_plasmids=("source_file", "nunique"),
    ).sort_values("total_hits", ascending=False)
    print(fam_counts.to_string())

    # By organism category
    print(f"\n--- IS hits by organism category ---")
    cat_counts = is_df.groupby(["organism_category", "is_family"]).size().reset_index(name="count")
    for cat in ["Gram-negative", "Gram-positive", "Acinetobacter", "Pseudomonas"]:
        sub = cat_counts[cat_counts["organism_category"] == cat].sort_values("count", ascending=False)
        if len(sub) > 0:
            print(f"\n  {cat}:")
            for _, row in sub.iterrows():
                print(f"    {row['is_family']}: {row['count']}")

    # By Inc/rep group
    print(f"\n--- IS hits by Inc/rep group ---")
    group_counts = is_df.groupby("inc_group").agg(
        total_IS=("source_file", "count"),
        plasmids_with_IS=("source_file", "nunique"),
    ).sort_values("total_IS", ascending=False)
    print(group_counts.to_string())

    # Detect composite transposons
    print(f"\n--- Composite Transposons ---")
    ct_df = detect_composite_transposons(is_df, AMR_FILE)
    if len(ct_df) > 0:
        ct_df.to_csv(os.path.join(OUT_DIR, "composite_transposons_curated.tsv"), sep="\t", index=False)
        print(f"Total composite transposons: {len(ct_df)}")

        # Top IS-AMR pairs
        is_amr = []
        for _, row in ct_df.iterrows():
            for gene in row["cargo_genes"].split(","):
                is_amr.append(f"{row['is_family']}-{gene}")
        pair_counts = Counter(is_amr).most_common(15)
        print("\nTop IS-AMR associations:")
        for pair, count in pair_counts:
            print(f"  {pair}: n={count}")

        # CT by group
        print("\nComposite transposons by Inc/rep group:")
        ct_by_group = ct_df.groupby("inc_group").size().sort_values(ascending=False)
        for g, n in ct_by_group.items():
            print(f"  {g}: {n}")
    else:
        print("No composite transposons detected.")

    # Save summary
    summary = {
        "total_plasmids": total_plasmids,
        "plasmids_with_IS": total_with_IS,
        "pct_with_IS": round(100 * total_with_IS / total_plasmids, 1),
        "total_IS_hits": len(is_df),
        "unique_IS_families": int(is_df["is_family"].nunique()),
        "top_IS_families": fam_counts.reset_index().to_dict(orient="records"),
        "composite_transposons": len(ct_df),
        "blast_params": "evalue=1e-20, perc_identity=85%, coverage>=70%, task=blastn",
        "reference_db": "20 curated IS elements from ISfinder/NCBI",
    }
    with open(os.path.join(OUT_DIR, "is_summary_curated.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nOutput saved to {OUT_DIR}")


if __name__ == "__main__":
    main()
