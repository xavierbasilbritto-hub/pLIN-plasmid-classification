#!/usr/bin/env python3
"""
Expanded pLIN outbreak validation — classify 57 new outbreak plasmids.

Uses the trained KNN classifier (data/inc_classifier.npz) and the full training
set pLIN assignments (output/pLIN_assignments.tsv) to:
1. Predict Inc type for each query plasmid
2. Find nearest-neighbour in the training set
3. Assign pLIN code from the nearest neighbour
4. Output results in the same format as outbreak_validation_pLIN_results.tsv
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

# ── Paths ─────────────────────────────────────────────────────────────────────

CLASSIFIER_PATH = os.path.join(BASE_DIR, "data", "inc_classifier.npz")
TRAINING_PLIN = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")
EXPANDED_DIR = os.path.join(os.path.dirname(__file__), "expanded_sequences")
DOWNLOAD_SUMMARY = os.path.join(EXPANDED_DIR, "download_summary.tsv")
OUTPUT_FILE = os.path.join(BASE_DIR, "output", "outbreak_validation_expanded_results.tsv")

# Metadata for each accession (from download script)
ACCESSION_META = {
    "GU595196": {"gene": "blaKPC-2", "country": "USA", "study": "Kitchel_2009_KPC_IncN"},
    "JN233704": {"gene": "blaKPC-2", "country": "USA", "study": "Chen_2012_KPC_IncFIA"},
    "CP004366": {"gene": "blaKPC-3", "country": "USA", "study": "Conlan_2014_NIH_KPC"},
    "CP004367": {"gene": "blaKPC-3", "country": "USA", "study": "Conlan_2014_NIH_KPC"},
    "CP019026": {"gene": "blaKPC-2", "country": "China", "study": "Li_2018_KPC_IncN"},
    "CP081510": {"gene": "blaKPC-2", "country": "Italy", "study": "Arcari_2023_KPC_outbreak"},
    "CP081509": {"gene": "blaKPC-2", "country": "Italy", "study": "Arcari_2023_KPC_outbreak"},
    "MH133192": {"gene": "blaKPC-2", "country": "Brazil", "study": "Andrade_2019_KPC_Brazil"},
    "MH234497": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234498": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234499": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234500": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234501": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234502": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234503": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234504": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234505": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234506": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234507": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234508": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234509": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "KX832926": {"gene": "blaNDM-1", "country": "Colombia", "study": "Rojas_2017_NDM_Colombia"},
    "KX832927": {"gene": "blaNDM-1", "country": "Colombia", "study": "Rojas_2017_NDM_Colombia"},
    "KX832928": {"gene": "blaNDM-1", "country": "Colombia", "study": "Rojas_2017_NDM_Colombia"},
    "KX832929": {"gene": "blaNDM-1", "country": "Colombia", "study": "Rojas_2017_NDM_Colombia"},
    "CP017672": {"gene": "blaNDM-1", "country": "Colombia", "study": "Rojas_2017_NDM_Colombia"},
    "JX104760": {"gene": "blaNDM-1", "country": "China", "study": "Ho_2012_NDM_China"},
    "MH985166": {"gene": "blaNDM-5", "country": "China", "study": "Li_2020_NDM5_China"},
    "MH985167": {"gene": "blaNDM-5", "country": "China", "study": "Li_2020_NDM5_China"},
    "MH985168": {"gene": "blaNDM-5", "country": "China", "study": "Li_2020_NDM5_China"},
    "MH985169": {"gene": "blaNDM-5", "country": "China", "study": "Li_2020_NDM5_China"},
    "MH985170": {"gene": "blaNDM-5", "country": "China", "study": "Li_2020_NDM5_China"},
    "MH985171": {"gene": "blaNDM-5", "country": "China", "study": "Li_2020_NDM5_China"},
    "JN626286": {"gene": "blaOXA-48", "country": "Turkey", "study": "Potron_2013_OXA48_Turkey"},
    "LR025097": {"gene": "blaOXA-48", "country": "Netherlands", "study": "Jousset_2019_OXA48_NL"},
    "LR025098": {"gene": "blaOXA-48", "country": "Netherlands", "study": "Jousset_2019_OXA48_NL"},
    "LR025100": {"gene": "blaOXA-48", "country": "Netherlands", "study": "Jousset_2019_OXA48_NL"},
    "LR025105": {"gene": "blaOXA-48", "country": "Netherlands", "study": "Jousset_2019_OXA48_NL"},
    "KP061858": {"gene": "blaOXA-48", "country": "France", "study": "Jousset_2019_OXA48_FR"},
    "MN783743": {"gene": "blaVIM-1", "country": "Italy", "study": "Arcari_2020_VIM_Italy"},
    "MN783744": {"gene": "blaVIM-1", "country": "Italy", "study": "Arcari_2020_VIM_Italy"},
    "MN783745": {"gene": "blaVIM-1", "country": "Italy", "study": "Arcari_2020_VIM_Italy"},
    "AB616660": {"gene": "blaIMP-6", "country": "Japan", "study": "Tada_2015_IMP_Japan"},
    "KP347127": {"gene": "mcr-1", "country": "China", "study": "Liu_2016_mcr1_discovery"},
    "KU761326": {"gene": "mcr-1", "country": "China", "study": "Zheng_2017_mcr1_China"},
    "KU761327": {"gene": "mcr-1", "country": "China", "study": "Zheng_2017_mcr1_China"},
    "KY075653": {"gene": "mcr-1", "country": "Europe", "study": "Hasman_2015_mcr1_Europe"},
    "KY075654": {"gene": "mcr-1", "country": "Europe", "study": "Hasman_2015_mcr1_Europe"},
    "CP016405": {"gene": "mcr-1", "country": "USA", "study": "McGann_2016_mcr1_USA"},
    "AY458016": {"gene": "blaCTX-M-15", "country": "India", "study": "Karim_2001_CTXM_India"},
    "EU935738": {"gene": "blaCTX-M-15", "country": "UK", "study": "Woodford_2009_CTXM_UK"},
    "EU935739": {"gene": "blaCTX-M-15", "country": "UK", "study": "Woodford_2009_CTXM_UK"},
    "EU935740": {"gene": "blaCTX-M-15", "country": "UK", "study": "Woodford_2009_CTXM_UK"},
    "CP009231": {"gene": "blaCTX-M-15", "country": "USA", "study": "Sheppard_2016_CTXM_USA"},
    "CP009232": {"gene": "blaCTX-M-15", "country": "USA", "study": "Sheppard_2016_CTXM_USA"},
    "CP009233": {"gene": "blaCTX-M-15", "country": "USA", "study": "Sheppard_2016_CTXM_USA"},
    "FN868832": {"gene": "blaCTX-M-15", "country": "France", "study": "Valverde_2009_CTXM_FR"},
}


def compute_kmer_vector(seq: str, k: int = 4) -> np.ndarray:
    """Compute normalised 4-mer frequency vector for a single sequence."""
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
    print("Expanded pLIN Outbreak Validation")
    print("=" * 70)

    # ── Load classifier ───────────────────────────────────────────────────
    print("\n[1/5] Loading KNN classifier ...")
    data = np.load(CLASSIFIER_PATH, allow_pickle=True)
    X_train = data["X"]          # (6998, 256)
    y_train = data["y"]          # (6998,) — integer labels
    group_names = data["group_names"]  # Inc type names
    print(f"  Training set: {X_train.shape[0]} plasmids, {len(group_names)} Inc types")

    # ── Load training pLIN assignments ────────────────────────────────────
    print("\n[2/5] Loading training pLIN assignments ...")
    df_plin = pd.read_csv(TRAINING_PLIN, sep="\t")
    print(f"  {len(df_plin)} training plasmids with pLIN codes")
    print(f"  {df_plin['pLIN'].nunique()} unique pLIN codes")

    # ── Load expanded sequences ───────────────────────────────────────────
    print("\n[3/5] Loading expanded outbreak sequences ...")
    fasta_files = sorted(glob.glob(os.path.join(EXPANDED_DIR, "*.fasta")))
    print(f"  Found {len(fasta_files)} FASTA files")

    queries = []
    for fpath in fasta_files:
        acc = os.path.basename(fpath).replace(".fasta", "")
        for rec in SeqIO.parse(fpath, "fasta"):
            seq_str = str(rec.seq)
            queries.append({
                "accession": acc,
                "header": rec.description[:120],
                "sequence": seq_str,
                "length_bp": len(seq_str),
            })
            break  # One sequence per file

    print(f"  Loaded {len(queries)} sequences")

    # ── Compute 4-mer vectors ─────────────────────────────────────────────
    print("\n[4/5] Computing 4-mer vectors and classifying ...")
    from scipy.spatial.distance import cosine as cosine_dist

    results = []
    for qi, q in enumerate(queries):
        acc = q["accession"]
        meta = ACCESSION_META.get(acc, {"gene": "unknown", "country": "unknown", "study": "unknown"})

        # Compute 4-mer vector
        q_vec = compute_kmer_vector(q["sequence"])

        # KNN classification — compute cosine distances to all training samples
        dists = np.array([cosine_dist(q_vec, X_train[j]) for j in range(X_train.shape[0])])

        # k=5 nearest neighbours
        nn_idx = np.argsort(dists)[:5]
        nn_dists = dists[nn_idx]
        nn_labels = y_train[nn_idx]

        # Distance-weighted voting
        weights = 1.0 / (nn_dists + 1e-10)
        vote_scores = {}
        for label, w in zip(nn_labels, weights):
            vote_scores[label] = vote_scores.get(label, 0.0) + w
        total_weight = sum(vote_scores.values())

        predicted_label = max(vote_scores, key=vote_scores.get)
        confidence = (vote_scores[predicted_label] / total_weight) * 100
        predicted_inc = group_names[predicted_label]

        # Nearest neighbour in training set
        nn_best_idx = nn_idx[0]
        nn_dist = nn_dists[0]
        nn_plasmid_id = df_plin.iloc[nn_best_idx]["plasmid_id"]
        nn_inc_type = df_plin.iloc[nn_best_idx]["inc_type"]
        nn_plin = df_plin.iloc[nn_best_idx]["pLIN"]

        # Secondary Inc type (if applicable)
        sorted_votes = sorted(vote_scores.items(), key=lambda x: -x[1])
        secondary_inc = ""
        if len(sorted_votes) > 1:
            sec_label, sec_weight = sorted_votes[1]
            sec_conf = (sec_weight / total_weight) * 100
            if sec_conf > 25:
                secondary_inc = f"{group_names[sec_label]} ({sec_conf:.1f}%)"

        results.append({
            "accession": acc,
            "study": meta["study"],
            "resistance_gene": meta["gene"],
            "country": meta["country"],
            "length_bp": q["length_bp"],
            "predicted_inc": str(predicted_inc),
            "confidence": round(confidence, 1),
            "secondary_inc": secondary_inc,
            "pLIN": nn_plin,
            "nn_plasmid": nn_plasmid_id,
            "nn_distance": round(float(nn_dist), 6),
            "nn_inc_type": nn_inc_type,
            "nn_plin": nn_plin,
        })

        status = "✓" if confidence >= 60 else "~"
        print(f"  [{qi+1}/{len(queries)}] {status} {acc}: {predicted_inc} ({confidence:.1f}%) "
              f"pLIN={nn_plin} d={nn_dist:.6f}")

    # ── Save results ──────────────────────────────────────────────────────
    print(f"\n[5/5] Saving results ...")
    df_results = pd.DataFrame(results)
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    df_results.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"  Saved to: {OUTPUT_FILE}")

    # ── Summary statistics ────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    print(f"  Total plasmids:         {len(df_results)}")
    print(f"  Countries:              {df_results['country'].nunique()}")
    print(f"  Studies:                {df_results['study'].nunique()}")
    print(f"  Resistance genes:       {df_results['resistance_gene'].nunique()}")
    print(f"  Unique pLIN L6 codes:   {df_results['pLIN'].nunique()}")

    high_conf = df_results[df_results["confidence"] >= 60]
    low_conf = df_results[df_results["confidence"] < 60]
    print(f"\n  High confidence (≥60%): {len(high_conf)} ({100*len(high_conf)/len(df_results):.1f}%)")
    print(f"  Low confidence (<60%):  {len(low_conf)} ({100*len(low_conf)/len(df_results):.1f}%)")

    print(f"\n  Mean NN distance:       {df_results['nn_distance'].mean():.6f}")
    print(f"  Median NN distance:     {df_results['nn_distance'].median():.6f}")
    print(f"  Min NN distance:        {df_results['nn_distance'].min():.6f}")
    print(f"  Max NN distance:        {df_results['nn_distance'].max():.6f}")

    # Per-gene summary
    print(f"\n  Per-resistance-gene breakdown:")
    for gene in sorted(df_results["resistance_gene"].unique()):
        sub = df_results[df_results["resistance_gene"] == gene]
        mean_conf = sub["confidence"].mean()
        mean_dist = sub["nn_distance"].mean()
        n_unique_plin = sub["pLIN"].nunique()
        print(f"    {gene:<16s}: n={len(sub):>2d}, mean_conf={mean_conf:.1f}%, "
              f"mean_nn_dist={mean_dist:.6f}, unique_pLIN={n_unique_plin}")

    # Outbreak clustering: plasmids sharing same pLIN within same study
    print(f"\n  Intra-study clustering (same pLIN within same study):")
    for study in sorted(df_results["study"].unique()):
        sub = df_results[df_results["study"] == study]
        if len(sub) < 2:
            continue
        plin_counts = sub["pLIN"].value_counts()
        shared = plin_counts[plin_counts > 1]
        if not shared.empty:
            for code, n in shared.items():
                print(f"    {study}: {n} plasmids share pLIN {code}")

    print("\n" + "=" * 70)
    print("Done!")


if __name__ == "__main__":
    main()
