#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Honest, leak-free validation of a contrastive-encoder feature transform
against the existing raw-4-mer KNN classifier for Inc/Rep group
classification.

Background: an internal discussion document (output/pLIN_Classifier_KNN_vs_
Encoder_Discussion.docx) reported a contrastive encoder beating KNN
(93.3% accuracy / 0.723 macro-F1 vs 91.1% / 0.666), but no training code,
saved model, or CV fold logs existed anywhere in this repository to verify
those numbers: they were hardcoded directly into a slide-deck generation
script with no underlying analysis to re-run. This script builds a real
encoder and runs a genuinely leak-free comparison, so the result is
independently reproducible from this file alone.

Leakage control (the discussion document's own open question #1): for
each of 5 stratified CV folds, the encoder is trained FROM SCRATCH on that
fold's training split only, and only ever sees the held-out test fold at
evaluation time, exactly mirroring how the existing KNN baseline's
cv_accuracy/cv_metrics in data/inc_classifier.npz were computed. Both
classifiers are evaluated identically: predict Inc/Rep group for each
held-out plasmid via distance-weighted k=5 nearest neighbour vote (cosine
distance): the only difference is which feature space the vote happens
in (raw 256-dim 4-mer frequencies for KNN, the encoder's learned
embedding for the encoder condition). This isolates the actual variable
of interest (does the learned embedding separate classes better than raw
composition) rather than confounding it with a different classification
rule.

Usage:
  python validate_encoder_vs_knn.py
  python validate_encoder_vs_knn.py --folds 5 --epochs 60 --embed-dim 64

Output:
  Prints per-fold and pooled accuracy/macro-F1/per-class F1 for both
  conditions to stdout.
  Writes output/encoder_vs_knn_validation_result.json (full methodology +
  results, machine-readable).
"""

import argparse
import json
import os
import time
from datetime import datetime, timezone

# torch must be imported before sklearn on this platform: importing
# sklearn first and torch second causes a native OpenMP/BLAS runtime
# conflict (both link their own copy) that segfaults on first tensor op
#: reproduced and confirmed during development of this script.
import torch
import torch.nn as nn
import torch.nn.functional as F

import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, f1_score
from sklearn.neighbors import KNeighborsClassifier

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CLASSIFIER_PATH = os.path.join(BASE_DIR, "data", "inc_classifier.npz")


class ContrastiveEncoder(nn.Module):
    """Small MLP embedding network: 256-dim 4-mer vector -> embed_dim
    L2-normalised embedding, trained with a supervised contrastive loss
    (same-class pairs pulled together, different-class pairs pushed
    apart): the architecture described in the discussion document
    (nonlinear embedding for cosine-KNN), reconstructed here since no
    saved implementation existed."""

    def __init__(self, in_dim=256, hidden_dim=128, embed_dim=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.ReLU(),
            nn.BatchNorm1d(hidden_dim),
            nn.Linear(hidden_dim, embed_dim),
        )

    def forward(self, x):
        z = self.net(x)
        return F.normalize(z, dim=1)


def supervised_contrastive_loss(embeddings, labels, temperature=0.1):
    """SupCon loss (Khosla et al. 2020): for each anchor, pull its
    embedding toward same-label embeddings and push away from
    different-label embeddings within the batch."""
    device = embeddings.device
    sim = torch.matmul(embeddings, embeddings.T) / temperature
    n = embeddings.shape[0]
    labels = labels.view(-1, 1)
    mask_pos = (labels == labels.T).float().to(device)
    mask_self = torch.eye(n, device=device)
    mask_pos = mask_pos - mask_self  # exclude self as its own positive

    logits_max, _ = sim.max(dim=1, keepdim=True)
    sim = sim - logits_max.detach()  # numerical stability
    exp_sim = torch.exp(sim) * (1 - mask_self)
    log_prob = sim - torch.log(exp_sim.sum(dim=1, keepdim=True) + 1e-12)

    n_pos = mask_pos.sum(dim=1)
    valid = n_pos > 0
    loss_per_anchor = -(mask_pos * log_prob).sum(dim=1)[valid] / n_pos[valid]
    return loss_per_anchor.mean()


def train_encoder(X_train, y_train, embed_dim=64, hidden_dim=128, epochs=60,
                   batch_size=256, lr=1e-3, device="cpu", seed=42):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = ContrastiveEncoder(in_dim=X_train.shape[1], hidden_dim=hidden_dim,
                                embed_dim=embed_dim).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    X_t = torch.tensor(X_train, dtype=torch.float32).to(device)
    y_t = torch.tensor(y_train, dtype=torch.long).to(device)
    n = X_t.shape[0]

    model.train()
    for epoch in range(epochs):
        perm = torch.randperm(n)
        total_loss = 0.0
        n_batches = 0
        for i in range(0, n, batch_size):
            idx = perm[i:i + batch_size]
            xb, yb = X_t[idx], y_t[idx]
            if len(torch.unique(yb)) < 2:
                continue  # SupCon needs at least 2 classes per batch to be meaningful
            emb = model(xb)
            loss = supervised_contrastive_loss(emb, yb)
            opt.zero_grad()
            loss.backward()
            opt.step()
            total_loss += loss.item()
            n_batches += 1
    model.eval()
    return model


def embed(model, X, device="cpu"):
    with torch.no_grad():
        X_t = torch.tensor(X, dtype=torch.float32).to(device)
        return model(X_t).cpu().numpy()


def evaluate_knn_vote(X_train, y_train, X_test, k=5, metric="cosine"):
    """Identical classification rule used for both conditions: k=5,
    cosine-distance, distance-weighted KNN vote: matching the production
    classifier's own configuration exactly."""
    clf = KNeighborsClassifier(n_neighbors=k, metric=metric, weights="distance")
    clf.fit(X_train, y_train)
    return clf.predict(X_test)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--embed-dim", type=int, default=64)
    ap.add_argument("--hidden-dim", type=int, default=128)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--output-json", default=os.path.join(BASE_DIR, "output", "encoder_vs_knn_validation_result.json"))
    args = ap.parse_args()

    print(f"Loading {CLASSIFIER_PATH} ...")
    data = np.load(CLASSIFIER_PATH, allow_pickle=True)
    X = data["X"].astype(np.float64)
    y = data["y"]
    group_names = [str(g) for g in data["group_names"]]
    print(f"  {X.shape[0]} samples, {X.shape[1]} features, {len(group_names)} classes")

    skf = StratifiedKFold(n_splits=args.folds, shuffle=True, random_state=args.seed)

    knn_preds_all, encoder_preds_all, y_true_all = [], [], []
    fold_results = []

    for fold_i, (train_idx, test_idx) in enumerate(skf.split(X, y)):
        t0 = time.time()
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Condition 1: raw 4-mer KNN (existing production classifier, no change)
        knn_pred = evaluate_knn_vote(X_train, y_train, X_test)

        # Condition 2: encoder trained FROM SCRATCH on this fold's train
        # split only (leak-free), then KNN vote in the resulting embedding
        # space. Encoder never sees X_test/y_test during training.
        model = train_encoder(X_train, y_train, embed_dim=args.embed_dim,
                               hidden_dim=args.hidden_dim, epochs=args.epochs,
                               seed=args.seed + fold_i)
        X_train_emb = embed(model, X_train)
        X_test_emb = embed(model, X_test)
        encoder_pred = evaluate_knn_vote(X_train_emb, y_train, X_test_emb)

        knn_acc = accuracy_score(y_test, knn_pred)
        knn_f1 = f1_score(y_test, knn_pred, average="macro", zero_division=0)
        enc_acc = accuracy_score(y_test, encoder_pred)
        enc_f1 = f1_score(y_test, encoder_pred, average="macro", zero_division=0)

        elapsed = time.time() - t0
        print(f"Fold {fold_i+1}/{args.folds} ({elapsed:.1f}s): "
              f"KNN acc={knn_acc:.4f} macroF1={knn_f1:.4f}  |  "
              f"Encoder acc={enc_acc:.4f} macroF1={enc_f1:.4f}")

        fold_results.append({
            "fold": fold_i, "n_test": len(test_idx),
            "knn_accuracy": knn_acc, "knn_macro_f1": knn_f1,
            "encoder_accuracy": enc_acc, "encoder_macro_f1": enc_f1,
        })
        knn_preds_all.append(knn_pred)
        encoder_preds_all.append(encoder_pred)
        y_true_all.append(y_test)

    y_true_all = np.concatenate(y_true_all)
    knn_preds_all = np.concatenate(knn_preds_all)
    encoder_preds_all = np.concatenate(encoder_preds_all)

    pooled_knn_acc = accuracy_score(y_true_all, knn_preds_all)
    pooled_knn_f1 = f1_score(y_true_all, knn_preds_all, average="macro", zero_division=0)
    pooled_enc_acc = accuracy_score(y_true_all, encoder_preds_all)
    pooled_enc_f1 = f1_score(y_true_all, encoder_preds_all, average="macro", zero_division=0)

    knn_per_class = f1_score(y_true_all, knn_preds_all, average=None, zero_division=0, labels=range(len(group_names)))
    enc_per_class = f1_score(y_true_all, encoder_preds_all, average=None, zero_division=0, labels=range(len(group_names)))

    print()
    print("=" * 70)
    print("POOLED CROSS-VALIDATED RESULTS (all folds combined)")
    print("=" * 70)
    print(f"  KNN (raw 4-mer):        accuracy={pooled_knn_acc:.4f}  macro-F1={pooled_knn_f1:.4f}")
    print(f"  Contrastive encoder:    accuracy={pooled_enc_acc:.4f}  macro-F1={pooled_enc_f1:.4f}")
    print(f"  Delta:                  accuracy={pooled_enc_acc - pooled_knn_acc:+.4f}  macro-F1={pooled_enc_f1 - pooled_knn_f1:+.4f}")
    print()
    print("Per-class F1 (KNN -> Encoder), classes with |delta| >= 0.05:")
    for i, name in enumerate(group_names):
        delta = enc_per_class[i] - knn_per_class[i]
        if abs(delta) >= 0.05:
            print(f"  {name:20s} n={int((y == i).sum()):5d}  KNN={knn_per_class[i]:.3f} -> Encoder={enc_per_class[i]:.3f}  ({delta:+.3f})")

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "methodology": (
            "Leak-free stratified k-fold CV. Both conditions use an "
            "identical k=5, cosine-distance, distance-weighted KNN vote "
            "as the classification rule; the only difference is the "
            "feature space (raw 256-dim 4-mer frequencies vs. a "
            "contrastive encoder's learned embedding). The encoder is "
            "trained from scratch on each fold's training split only and "
            "never sees the held-out test fold during training."
        ),
        "n_folds": args.folds, "encoder_epochs": args.epochs,
        "encoder_embed_dim": args.embed_dim, "encoder_hidden_dim": args.hidden_dim,
        "seed": args.seed,
        "pooled_results": {
            "knn_accuracy": pooled_knn_acc, "knn_macro_f1": pooled_knn_f1,
            "encoder_accuracy": pooled_enc_acc, "encoder_macro_f1": pooled_enc_f1,
            "delta_accuracy": pooled_enc_acc - pooled_knn_acc,
            "delta_macro_f1": pooled_enc_f1 - pooled_knn_f1,
        },
        "per_fold_results": fold_results,
        "per_class_f1": {
            group_names[i]: {"knn": float(knn_per_class[i]), "encoder": float(enc_per_class[i])}
            for i in range(len(group_names))
        },
        "comparison_to_previously_unverified_claim": {
            "unverified_claim_source": "output/pLIN_Classifier_KNN_vs_Encoder_Discussion.docx (hardcoded in create_classifier_comparison_pptx.py, no training code or CV logs existed)",
            "unverified_claim_accuracy": 0.933, "unverified_claim_macro_f1": 0.723,
            "this_scripts_genuine_result_accuracy": pooled_enc_acc,
            "this_scripts_genuine_result_macro_f1": pooled_enc_f1,
        },
    }
    os.makedirs(os.path.dirname(args.output_json), exist_ok=True)
    with open(args.output_json, "w") as f:
        json.dump(output, f, indent=2)
    print()
    print(f"Full results written to {args.output_json}")


if __name__ == "__main__":
    main()
