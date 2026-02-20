#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Train an Inc-group classifier from training data using k-mer profiles.

Approach: Compute 4-mer frequency vectors for all training plasmids,
then train a K-Nearest Neighbors classifier (cosine distance, k=5).
Also saves centroid profiles for reference.

Usage:  python build_inc_centroids.py
Output: data/inc_classifier.npz  (vectors + labels + centroids)
"""

import os
import sys
import glob
import numpy as np
from itertools import product as iter_product
from Bio import SeqIO

TRAINING_DIR = os.path.join(os.path.dirname(__file__), "plasmid_sequences_for_training")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "data")

K = 4
BASES = "ACGT"
KMERS = ["".join(p) for p in iter_product(BASES, repeat=K)]
KMER_INDEX = {km: i for i, km in enumerate(KMERS)}
N_KMERS = len(KMERS)  # 256


def kmer_vector(sequence):
    """Compute normalised 4-mer frequency vector for a single sequence."""
    seq = sequence.upper()
    counts = np.zeros(N_KMERS, dtype=np.float64)
    for i in range(len(seq) - K + 1):
        kmer = seq[i:i + K]
        if kmer in KMER_INDEX:
            counts[KMER_INDEX[kmer]] += 1
    total = counts.sum()
    if total > 0:
        counts /= total
    return counts


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    inc_groups = {}
    for inc_dir in sorted(glob.glob(os.path.join(TRAINING_DIR, "*", "fastas"))):
        inc_name = os.path.basename(os.path.dirname(inc_dir))
        fasta_files = sorted(glob.glob(os.path.join(inc_dir, "*.fasta")))
        if not fasta_files:
            continue
        inc_groups[inc_name] = fasta_files

    if not inc_groups:
        print("ERROR: No training data found in", TRAINING_DIR)
        sys.exit(1)

    all_vectors = []
    all_labels = []
    group_names = sorted(inc_groups.keys())
    label_map = {name: i for i, name in enumerate(group_names)}

    for inc_name in group_names:
        fasta_files = inc_groups[inc_name]
        print(f"\nProcessing {inc_name} ({len(fasta_files)} plasmids)...")
        count = 0
        for i, fpath in enumerate(fasta_files):
            if (i + 1) % 500 == 0 or i == 0:
                print(f"  [{i+1}/{len(fasta_files)}]")
            try:
                for rec in SeqIO.parse(fpath, "fasta"):
                    vec = kmer_vector(str(rec.seq))
                    all_vectors.append(vec)
                    all_labels.append(label_map[inc_name])
                    count += 1
                    break  # one sequence per file
            except Exception as e:
                print(f"  WARNING: skipping {fpath}: {e}")
        print(f"  → {count} vectors")

    X = np.array(all_vectors, dtype=np.float32)
    y = np.array(all_labels, dtype=np.int32)
    print(f"\nTotal training data: {X.shape[0]} sequences, {X.shape[1]} features")

    # Compute centroids per group
    centroids = np.zeros((len(group_names), N_KMERS), dtype=np.float64)
    for i, name in enumerate(group_names):
        mask = (y == i)
        centroids[i] = X[mask].mean(axis=0)

    # Cross-validation accuracy check
    import json
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.metrics import (classification_report, confusion_matrix,
                                 accuracy_score, precision_recall_fscore_support)

    # Dynamic k: ensure k < smallest class size
    min_class_size = min(np.bincount(y))
    k = min(5, max(1, min_class_size - 1))
    print(f"\nSmallest class has {min_class_size} samples → using k={k}")

    knn = KNeighborsClassifier(n_neighbors=k, metric="cosine", weights="distance")

    # Stratified CV for balanced evaluation across all Inc groups
    cv_metrics_json = "{}"
    conf_matrix = np.zeros((len(group_names), len(group_names)), dtype=np.int32)
    cv_accuracy = 0.0
    confusion_pairs_list = []

    n_splits = min(5, min_class_size)
    if n_splits < 2:
        print("WARNING: Some classes too small for cross-validation, skipping CV")
    else:
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
        y_pred = cross_val_predict(knn, X, y, cv=cv)
        print(f"\n{n_splits}-fold Stratified CV — Per-class metrics:")
        print(classification_report(y, y_pred, target_names=group_names, digits=4))

        # Compute per-class metrics and save
        cv_accuracy = float(accuracy_score(y, y_pred))
        prec, rec, f1, support = precision_recall_fscore_support(y, y_pred, zero_division=0)
        per_class = {}
        for i, name in enumerate(group_names):
            per_class[name] = {
                "precision": round(float(prec[i]), 4),
                "recall": round(float(rec[i]), 4),
                "f1": round(float(f1[i]), 4),
                "support": int(support[i]),
            }
        cv_metrics_json = json.dumps(per_class)

        # Confusion matrix
        conf_matrix = confusion_matrix(y, y_pred)

        # Identify confusion pairs (off-diagonal > 5% of either class)
        for i in range(len(group_names)):
            for j in range(i + 1, len(group_names)):
                err_ij = conf_matrix[i, j]
                err_ji = conf_matrix[j, i]
                total_err = err_ij + err_ji
                pct_i = 100.0 * total_err / max(support[i], 1)
                pct_j = 100.0 * total_err / max(support[j], 1)
                if pct_i > 5 or pct_j > 5:
                    confusion_pairs_list.append(
                        f"{group_names[i]}|{group_names[j]}|{total_err}|{pct_i:.1f}|{pct_j:.1f}")

        print(f"\nOverall CV accuracy: {cv_accuracy:.4f}")
        if confusion_pairs_list:
            print(f"Confusion pairs (>5% error): {len(confusion_pairs_list)}")
            for cp in confusion_pairs_list:
                parts = cp.split("|")
                print(f"  {parts[0]} <-> {parts[1]}: {parts[2]} errors ({parts[3]}%/{parts[4]}%)")

    print("\nTraining set composition:")
    for i, name in enumerate(group_names):
        mask = y == i
        print(f"  {name}: {mask.sum()} samples")

    # Save everything (training vectors + labels + metadata + CV metrics)
    out_path = os.path.join(OUTPUT_DIR, "inc_classifier.npz")
    np.savez_compressed(
        out_path,
        X=X,
        y=y,
        group_names=np.array(group_names),
        centroids=centroids,
        kmers=np.array(KMERS),
        cv_metrics=np.array([cv_metrics_json]),
        confusion_matrix=conf_matrix,
        cv_accuracy=np.array([cv_accuracy]),
        confusion_pairs=np.array(confusion_pairs_list if confusion_pairs_list else [""]),
    )
    size_mb = os.path.getsize(out_path) / 1024 / 1024
    print(f"\nSaved: {out_path} ({size_mb:.1f} MB)")
    print(f"Inc groups: {group_names}")

    # Also keep the old centroid file for backward compat
    old_path = os.path.join(OUTPUT_DIR, "inc_centroids.npz")
    np.savez_compressed(
        old_path,
        group_names=np.array(group_names),
        centroids=centroids,
        stdevs=np.array([X[y == i].std(axis=0) for i in range(len(group_names))]),
        kmers=np.array(KMERS),
    )
    print(f"Also updated: {old_path}")


if __name__ == "__main__":
    main()
