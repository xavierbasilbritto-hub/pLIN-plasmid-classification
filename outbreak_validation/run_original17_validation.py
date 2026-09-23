#!/usr/bin/env python3
"""
Original 17-plasmid pLIN outbreak validation — reconstructed generator.

The script that originally produced output/outbreak_validation_pLIN_results.tsv
was never committed to the repository; only its (now-stale, pre-pipeline-fix)
output survived. This script reconstructs it using the identical
nearest-neighbour query methodology already used and verified in
run_expanded_validation.py (same KNN classifier, same k=5 distance-weighted
voting, same cosine-distance nearest-neighbour pLIN assignment), applied to
the 17 original outbreak FASTA files in outbreak_validation/*.fasta, against
the corrected, deterministic output/pLIN_assignments.tsv.

Metadata (study, plasmid_name, expected_inc, resistance_gene) is carried
over unchanged from the original stale results file — only the
classification (predicted_inc, confidence, pLIN, nn_plasmid, nn_distance)
is recomputed against the corrected training database.
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
from itertools import product as iter_product
from Bio import SeqIO

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BASE_DIR)

CLASSIFIER_PATH = os.path.join(BASE_DIR, "data", "inc_classifier.npz")
TRAINING_PLIN = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")
ORIGINAL_DIR = os.path.dirname(__file__)  # the 17 FASTAs sit directly here
OUTPUT_FILE = os.path.join(BASE_DIR, "output", "outbreak_validation_pLIN_results.tsv")

# Metadata carried over verbatim from the original (stale) results file —
# these are ground-truth facts about the source studies/plasmids, not
# classifier output, so they are unaffected by the pLIN renumbering fix.
ACCESSION_META = {
    "CP022533": {"study": "Study6_IMP4_IncHI2_Australia", "plasmid_name": "pMS7884A",
                 "expected_inc": "IncHI2", "resistance_gene": "blaIMP-4"},
    "CP104940": {"study": "Study1_KPC2_IncN_Germany", "plasmid_name": "pKV30046-KPC2",
                 "expected_inc": "IncN", "resistance_gene": "blaKPC-2"},
    "CP104944": {"study": "Study1_KPC2_IncN_Germany", "plasmid_name": "pKP37361-KPC2",
                 "expected_inc": "IncN", "resistance_gene": "blaKPC-2"},
    "CP104949": {"study": "Study1_KPC2_IncN_Germany", "plasmid_name": "pEClo_Surv151-KPC2",
                 "expected_inc": "IncN", "resistance_gene": "blaKPC-2"},
    "MN542377": {"study": "Study4_KPC2_Singapore", "plasmid_name": "pKPC2_sg1",
                 "expected_inc": "IncP/other", "resistance_gene": "blaKPC-2"},
    "MN657241": {"study": "Study3_NDM1_Germany", "plasmid_name": "pCF104a-T3",
                 "expected_inc": "IncAC2/IncN", "resistance_gene": "blaNDM-1"},
    "MN657242": {"study": "Study3_NDM1_Germany", "plasmid_name": "pEC405a-T3",
                 "expected_inc": "IncAC2/IncN", "resistance_gene": "blaNDM-1"},
    "MN657243": {"study": "Study3_NDM1_Germany", "plasmid_name": "pEC744-T5",
                 "expected_inc": "IncAC2/IncN", "resistance_gene": "blaNDM-1"},
    "MN657244": {"study": "Study3_NDM1_Germany", "plasmid_name": "pEC6332-T3",
                 "expected_inc": "IncAC2/IncN", "resistance_gene": "blaNDM-1"},
    "MN657245": {"study": "Study3_NDM1_Germany", "plasmid_name": "pEC6332-T6",
                 "expected_inc": "IncAC2/IncN", "resistance_gene": "blaNDM-1"},
    "MN657246": {"study": "Study3_NDM1_Germany", "plasmid_name": "pEC6332-T7",
                 "expected_inc": "IncAC2/IncN", "resistance_gene": "blaNDM-1"},
    "MN657247": {"study": "Study3_NDM1_Germany", "plasmid_name": "pECl-T3",
                 "expected_inc": "IncAC2/IncN", "resistance_gene": "blaNDM-1"},
    "MN657248": {"study": "Study3_NDM1_Germany", "plasmid_name": "pKP15-T2",
                 "expected_inc": "IncAC2/IncN", "resistance_gene": "blaNDM-1"},
    "MN657249": {"study": "Study3_NDM1_Germany", "plasmid_name": "pKP39-T3",
                 "expected_inc": "IncAC2/IncN", "resistance_gene": "blaNDM-1"},
    "MN657250": {"study": "Study3_NDM1_Germany", "plasmid_name": "pKP39-T4",
                 "expected_inc": "IncAC2/IncN", "resistance_gene": "blaNDM-1"},
    "MN657251": {"study": "Study3_NDM1_Germany", "plasmid_name": "pKPC-2",
                 "expected_inc": "IncN", "resistance_gene": "blaKPC-2"},
    "MN657252": {"study": "Study3_NDM1_Germany", "plasmid_name": "pPS-T1",
                 "expected_inc": "IncAC2/IncN", "resistance_gene": "blaNDM-1"},
}


def compute_kmer_vector(seq: str, k: int = 4) -> np.ndarray:
    """Compute normalised 4-mer frequency vector — identical to
    run_expanded_validation.py::compute_kmer_vector for consistency."""
    bases = "ACGT"
    all_kmers = ["".join(p) for p in iter_product(bases, repeat=k)]
    vec = np.zeros(len(all_kmers), dtype=np.float64)
    seq = seq.upper()
    total = max(len(seq) - k + 1, 1)
    for ki, kmer in enumerate(all_kmers):
        vec[ki] = seq.count(kmer) / total
    return vec


def main():
    print("=" * 70)
    print("Original 17-Plasmid pLIN Outbreak Validation (reconstructed)")
    print("=" * 70)

    print("\n[1/5] Loading KNN classifier ...")
    data = np.load(CLASSIFIER_PATH, allow_pickle=True)
    X_train = data["X"]
    y_train = data["y"]
    group_names = data["group_names"]
    print(f"  Training set: {X_train.shape[0]} plasmids, {len(group_names)} Inc types")

    print("\n[2/5] Loading training pLIN assignments ...")
    df_plin = pd.read_csv(TRAINING_PLIN, sep="\t")
    print(f"  {len(df_plin)} training plasmids with pLIN codes")
    print(f"  {df_plin['pLIN'].nunique()} unique pLIN codes")

    print("\n[3/5] Loading original 17 outbreak sequences ...")
    fasta_files = sorted(
        f for f in glob.glob(os.path.join(ORIGINAL_DIR, "*.fasta"))
        if os.path.splitext(os.path.basename(f))[0] in ACCESSION_META
    )
    print(f"  Found {len(fasta_files)} FASTA files")
    if len(fasta_files) != len(ACCESSION_META):
        missing = set(ACCESSION_META) - {os.path.splitext(os.path.basename(f))[0] for f in fasta_files}
        print(f"  WARNING: missing FASTA files for: {sorted(missing)}")

    queries = []
    for fpath in fasta_files:
        acc = os.path.splitext(os.path.basename(fpath))[0]
        for rec in SeqIO.parse(fpath, "fasta"):
            seq_str = str(rec.seq)
            queries.append({
                "accession": acc,
                "sequence": seq_str,
                "length_bp": len(seq_str),
            })
            break  # One sequence per file

    print(f"  Loaded {len(queries)} sequences")

    print("\n[4/5] Computing 4-mer vectors and classifying ...")
    from scipy.spatial.distance import cosine as cosine_dist

    results = []
    for qi, q in enumerate(queries):
        acc = q["accession"]
        meta = ACCESSION_META[acc]

        q_vec = compute_kmer_vector(q["sequence"])

        dists = np.array([cosine_dist(q_vec, X_train[j]) for j in range(X_train.shape[0])])

        nn_idx = np.argsort(dists)[:5]
        nn_dists = dists[nn_idx]
        nn_labels = y_train[nn_idx]

        weights = 1.0 / (nn_dists + 1e-10)
        vote_scores = {}
        for label, w in zip(nn_labels, weights):
            vote_scores[label] = vote_scores.get(label, 0.0) + w
        total_weight = sum(vote_scores.values())

        predicted_label = max(vote_scores, key=vote_scores.get)
        confidence = (vote_scores[predicted_label] / total_weight) * 100
        predicted_inc = group_names[predicted_label]

        nn_best_idx = nn_idx[0]
        nn_dist = nn_dists[0]
        nn_plasmid_id = df_plin.iloc[nn_best_idx]["plasmid_id"]
        nn_inc_type = df_plin.iloc[nn_best_idx]["inc_type"]
        nn_plin = df_plin.iloc[nn_best_idx]["pLIN"]

        results.append({
            "accession": acc,
            "study": meta["study"],
            "plasmid_name": meta["plasmid_name"],
            "length_bp": q["length_bp"],
            "predicted_inc": str(predicted_inc),
            "confidence": round(confidence, 1),
            "expected_inc": meta["expected_inc"],
            "pLIN": nn_plin,
            "nn_plasmid": nn_plasmid_id,
            "nn_distance": round(float(nn_dist), 6),
            "nn_inc_type": nn_inc_type,
            "nn_plin": nn_plin,
            "resistance_gene": meta["resistance_gene"],
        })

        status = "✓" if confidence >= 60 else "~"
        print(f"  [{qi+1}/{len(queries)}] {status} {acc}: {predicted_inc} ({confidence:.1f}%) "
              f"pLIN={nn_plin} d={nn_dist:.6f}")

    print(f"\n[5/5] Saving results ...")
    df_results = pd.DataFrame(results)
    # Preserve original column order
    df_results = df_results[["accession", "study", "plasmid_name", "length_bp",
                              "predicted_inc", "confidence", "expected_inc", "pLIN",
                              "nn_plasmid", "nn_distance", "nn_inc_type", "nn_plin",
                              "resistance_gene"]]
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    df_results.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"  Saved to: {OUTPUT_FILE}")

    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    print(f"  Total plasmids: {len(df_results)}")
    print(f"  Studies: {df_results['study'].nunique()}")
    print(f"  Unique pLIN codes: {df_results['pLIN'].nunique()}")
    match = (df_results["predicted_inc"].str.split("/").str[0] ==
              df_results["expected_inc"].str.split("/").str[0])
    print(f"  predicted_inc matches expected_inc (loosely): {match.sum()}/{len(df_results)}")

    print("\nDone!")


if __name__ == "__main__":
    main()
