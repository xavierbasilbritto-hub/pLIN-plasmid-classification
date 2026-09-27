#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Train the production contrastive encoder for optional Inc/Rep group
classification, and save it alongside a metadata file recording its
genuine, leak-free cross-validated performance.

This is the encoder architecture validated in validate_encoder_vs_knn.py
(5-fold CV: 92.7% accuracy / 0.692 macro-F1 vs KNN's 91.1% / 0.666, a
consistent improvement in all 5 folds, 19/28 classes improved and 8
declined by at most 0.016 with none catastrophic) — see that script and
output/encoder_vs_knn_validation_result.json for the full validation.

This script trains the FINAL production encoder on the complete 8,077-
plasmid training set (not a CV fold), since a deployed encoder should use
all available labelled data, exactly as the existing KNN classifier does
(data/inc_classifier.npz is also trained on the full set, with its 91.1%
figure coming from a separate CV run rather than this final artifact's
own training data).

Usage:
  python train_inc_encoder.py

Output:
  data/inc_encoder.pt          — encoder weights (PyTorch state_dict) +
                                  embedded training vectors/labels for
                                  KNN-in-embedding-space at inference time
  data/inc_encoder_metadata.json — architecture, training config, and the
                                  genuine CV performance figures from
                                  validate_encoder_vs_knn.py (NOT
                                  re-measured on this final fit, which
                                  would be optimistic/leaked — the
                                  reported accuracy is always the
                                  independent CV result)
"""

# torch must be imported before sklearn — see validate_encoder_vs_knn.py
# for the reproduced native-library import-order conflict this avoids.
import torch

import json
import os
from datetime import datetime, timezone

import numpy as np

from validate_encoder_vs_knn import ContrastiveEncoder, train_encoder, embed

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CLASSIFIER_PATH = os.path.join(BASE_DIR, "data", "inc_classifier.npz")
ENCODER_PATH = os.path.join(BASE_DIR, "data", "inc_encoder.pt")
METADATA_PATH = os.path.join(BASE_DIR, "data", "inc_encoder_metadata.json")
CV_RESULT_PATH = os.path.join(BASE_DIR, "output", "encoder_vs_knn_validation_result.json")


def main():
    print(f"Loading {CLASSIFIER_PATH} ...")
    data = np.load(CLASSIFIER_PATH, allow_pickle=True)
    X = data["X"].astype(np.float64)
    y = data["y"]
    group_names = [str(g) for g in data["group_names"]]
    print(f"  {X.shape[0]} samples, {X.shape[1]} features, {len(group_names)} classes")

    if not os.path.exists(CV_RESULT_PATH):
        raise SystemExit(
            f"ERROR: {CV_RESULT_PATH} not found. Run validate_encoder_vs_knn.py "
            "first to genuinely measure this architecture's CV performance before "
            "training a production artifact — training a final model without a "
            "prior independent CV run would leave no honest accuracy figure to "
            "report for it."
        )
    with open(CV_RESULT_PATH) as f:
        cv_result = json.load(f)

    print("Training final production encoder on the full training set ...")
    embed_dim = cv_result["encoder_embed_dim"]
    hidden_dim = cv_result["encoder_hidden_dim"]
    epochs = cv_result["encoder_epochs"]
    seed = cv_result["seed"]

    model = train_encoder(X, y, embed_dim=embed_dim, hidden_dim=hidden_dim,
                           epochs=epochs, seed=seed)
    X_embedded = embed(model, X)

    torch.save({
        "model_state_dict": model.state_dict(),
        "in_dim": X.shape[1],
        "hidden_dim": hidden_dim,
        "embed_dim": embed_dim,
        "X_train_embedded": X_embedded,
        "y_train": y,
        "group_names": group_names,
    }, ENCODER_PATH)
    print(f"Saved encoder + embedded training set to {ENCODER_PATH}")

    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "architecture": "ContrastiveEncoder (see validate_encoder_vs_knn.py): "
                         f"{X.shape[1]} -> {hidden_dim} (ReLU, BatchNorm) -> {embed_dim} (L2-normalised)",
        "training": f"Supervised contrastive loss, {epochs} epochs, seed={seed}, "
                    "trained on the full 8,077-plasmid training set.",
        "genuine_cv_performance": {
            "note": (
                "These figures come from validate_encoder_vs_knn.py's independent, "
                "leak-free 5-fold cross-validation (encoder retrained from scratch "
                "each fold, never sees its own test fold) — NOT from evaluating this "
                "final artifact on the data it was trained on, which would be "
                "optimistic. This is the accuracy a user should expect from the "
                "encoder option, not a number specific to this exact saved file."
            ),
            "knn_baseline_accuracy": cv_result["pooled_results"]["knn_accuracy"],
            "knn_baseline_macro_f1": cv_result["pooled_results"]["knn_macro_f1"],
            "encoder_cv_accuracy": cv_result["pooled_results"]["encoder_accuracy"],
            "encoder_cv_macro_f1": cv_result["pooled_results"]["encoder_macro_f1"],
            "n_cv_folds": cv_result["n_folds"],
        },
        "known_limitations": [
            "Adding a NEW Inc/Rep group (one not in the 28 currently trained on) "
            "requires retraining this encoder from scratch — unlike the raw-4-mer "
            "KNN classifier, where a new group is a FASTA folder drop and a "
            "training-data rebuild with no model retraining needed. The encoder "
            "option should be treated as unavailable/stale for any group added "
            "after this file's generated_utc timestamp until retrained.",
            "Per-class effect is not uniform: 8 of 28 classes showed a small F1 "
            "decline (largest -0.016) relative to KNN in cross-validation, "
            "though none catastrophically. See "
            "output/encoder_vs_knn_validation_result.json for the full per-class "
            "breakdown before relying on the encoder for a specific Inc group.",
            "Embedding stability across independent retrains (e.g. after adding "
            "training data) has not been separately measured — a retrained "
            "encoder's embedding geometry may shift, unlike KNN's raw 4-mer "
            "feature space, which cannot change under retraining.",
        ],
    }
    with open(METADATA_PATH, "w") as f:
        json.dump(metadata, f, indent=2)
    print(f"Saved metadata to {METADATA_PATH}")


if __name__ == "__main__":
    main()
