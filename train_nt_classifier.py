#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Train Inc-group and AMR classifiers using Nucleotide Transformer embeddings.

This script extracts embeddings from the pre-trained Nucleotide Transformer
genomic language model and trains lightweight linear probes for:
  1. Inc group classification (from labeled training sequences)
  2. AMR drug-class prediction (from AMRFinderPlus results, if available)

The trained probes are saved as pickle files for use in the pLIN Streamlit app.

Prerequisites:
    pip install transformers torch joblib

Usage:
    python train_nt_classifier.py [--model MODEL_NAME] [--device cpu|cuda|mps]

Output:
    data/nt_inc_probe.pkl       — Inc group classifier (LogisticRegression)
    data/nt_amr_probe.pkl       — AMR class predictor (optional)
    data/nt_embeddings.npz      — Cached embeddings for reuse
"""

import os
import sys
import glob
import argparse
import numpy as np
import joblib
from Bio import SeqIO

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
TRAINING_DIR = os.path.join(SCRIPT_DIR, "plasmid_sequences_for_training")
CLASSIFIER_PATH = os.path.join(DATA_DIR, "inc_classifier.npz")
AMR_RESULTS_PATH = os.path.join(SCRIPT_DIR, "output", "integrated",
                                 "pLIN_AMR_integrated.tsv")

NT_MODELS = {
    "50m": "InstaDeepAI/nucleotide-transformer-v2-50m-multi-species",
    "100m": "InstaDeepAI/nucleotide-transformer-v2-100m-multi-species",
    "250m": "InstaDeepAI/nucleotide-transformer-v2-250m-multi-species",
    "500m": "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
}

CHUNK_SIZE = 5000   # bp per chunk
STRIDE = 2500       # overlap stride


def extract_nt_embedding(sequence, tokenizer, model, device, chunk_size=CHUNK_SIZE,
                         stride=STRIDE):
    """Extract mean-pooled Nucleotide Transformer embedding for a sequence.

    Long sequences are split into overlapping chunks, embedded independently,
    and the chunk embeddings are averaged to produce a single vector.
    """
    import torch

    seq = sequence.upper()
    chunks = []
    for start in range(0, len(seq), stride):
        chunk = seq[start:start + chunk_size]
        if len(chunk) < 100:
            continue
        chunks.append(chunk)
    if not chunks:
        chunks = [seq[:chunk_size] if len(seq) >= 100 else seq]

    embeddings = []
    for chunk in chunks:
        tokens = tokenizer(chunk, return_tensors="pt", padding=True,
                           truncation=True, max_length=1000)
        tokens = {k: v.to(device) for k, v in tokens.items()}

        with torch.no_grad():
            outputs = model(**tokens, output_hidden_states=True)

        hidden = outputs.hidden_states[-1]  # last layer
        mask = tokens["attention_mask"].unsqueeze(-1).float()
        mean_emb = (hidden * mask).sum(dim=1) / mask.sum(dim=1)
        embeddings.append(mean_emb.squeeze().cpu().numpy())

    return np.mean(embeddings, axis=0).astype(np.float32)


def load_training_sequences():
    """Load training plasmid sequences and Inc group labels from the training directory."""
    inc_groups = {}
    for inc_dir in sorted(glob.glob(os.path.join(TRAINING_DIR, "*", "fastas"))):
        inc_name = os.path.basename(os.path.dirname(inc_dir))
        fasta_files = sorted(glob.glob(os.path.join(inc_dir, "*.fasta")))
        if fasta_files:
            inc_groups[inc_name] = fasta_files

    if not inc_groups:
        # Fallback: use inc_classifier.npz metadata if no FASTA files available
        print("No training FASTA files found. Cannot extract NT embeddings.")
        print(f"Expected directory: {TRAINING_DIR}/{{IncGroup}}/fastas/*.fasta")
        return None, None, None

    sequences = []
    labels = []
    filenames = []
    group_names = sorted(inc_groups.keys())
    label_map = {name: i for i, name in enumerate(group_names)}

    for inc_name in group_names:
        fasta_files = inc_groups[inc_name]
        print(f"Loading {inc_name}: {len(fasta_files)} files...")
        for fpath in fasta_files:
            try:
                for rec in SeqIO.parse(fpath, "fasta"):
                    sequences.append(str(rec.seq))
                    labels.append(label_map[inc_name])
                    filenames.append(os.path.basename(fpath))
                    break  # one sequence per file
            except Exception as e:
                print(f"  WARNING: skipping {fpath}: {e}")

    print(f"Loaded {len(sequences)} sequences from {len(group_names)} Inc groups")
    return sequences, np.array(labels), group_names


def load_amr_labels(filenames):
    """Load AMR class labels from integrated results (if available)."""
    import pandas as pd

    if not os.path.exists(AMR_RESULTS_PATH):
        print(f"No AMR results found at {AMR_RESULTS_PATH} — skipping AMR probe training")
        return None, None

    df = pd.read_csv(AMR_RESULTS_PATH, sep="\t")
    if "AMR_classes" not in df.columns or "source_file" not in df.columns:
        print("AMR results missing required columns — skipping AMR probe training")
        return None, None

    # Build multi-label AMR matrix
    all_classes = set()
    for classes_str in df["AMR_classes"].dropna():
        for cls in classes_str.split(";"):
            cls = cls.strip()
            if cls:
                all_classes.add(cls)

    if not all_classes:
        return None, None

    class_list = sorted(all_classes)
    class_to_idx = {c: i for i, c in enumerate(class_list)}

    # Map filenames to AMR labels
    file_to_amr = {}
    for _, row in df.iterrows():
        sf = str(row.get("source_file", ""))
        classes_str = str(row.get("AMR_classes", ""))
        amr_vec = np.zeros(len(class_list), dtype=np.int32)
        for cls in classes_str.split(";"):
            cls = cls.strip()
            if cls in class_to_idx:
                amr_vec[class_to_idx[cls]] = 1
        # Match on basename without extension
        for ext in [".fasta", ".fa", ".fna"]:
            sf_clean = sf.replace(ext, "")
            file_to_amr[sf_clean] = amr_vec
            file_to_amr[sf] = amr_vec

    amr_labels = []
    matched = 0
    for fn in filenames:
        fn_clean = fn.replace(".fasta", "").replace(".fa", "").replace(".fna", "")
        if fn_clean in file_to_amr:
            amr_labels.append(file_to_amr[fn_clean])
            matched += 1
        elif fn in file_to_amr:
            amr_labels.append(file_to_amr[fn])
            matched += 1
        else:
            amr_labels.append(np.zeros(len(class_list), dtype=np.int32))

    print(f"Matched {matched}/{len(filenames)} sequences to AMR labels "
          f"({len(class_list)} drug classes)")

    if matched < 10:
        print("Too few AMR matches — skipping AMR probe training")
        return None, None

    return np.array(amr_labels), class_list


def main():
    parser = argparse.ArgumentParser(
        description="Train NT-based Inc group and AMR classifiers")
    parser.add_argument("--model", default="50m",
                        choices=list(NT_MODELS.keys()),
                        help="NT model size (default: 50m)")
    parser.add_argument("--device", default=None,
                        choices=["cpu", "cuda", "mps"],
                        help="Compute device (auto-detected if not specified)")
    parser.add_argument("--batch-size", type=int, default=1,
                        help="Sequences per batch for embedding extraction")
    parser.add_argument("--max-sequences", type=int, default=None,
                        help="Limit number of sequences (for testing)")
    args = parser.parse_args()

    # ── Device setup ──
    import torch
    if args.device:
        device = torch.device(args.device)
    elif torch.cuda.is_available():
        device = torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        device = torch.device("mps")
    else:
        device = torch.device("cpu")
    print(f"Device: {device}")

    # ── Load model ──
    model_name = NT_MODELS[args.model]
    print(f"\nLoading Nucleotide Transformer: {model_name}")

    from transformers import AutoTokenizer, AutoModelForMaskedLM

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForMaskedLM.from_pretrained(model_name)
    model.eval()
    model.to(device)
    print(f"Model loaded ({sum(p.numel() for p in model.parameters()) / 1e6:.0f}M parameters)")

    hidden_size = model.config.hidden_size
    print(f"Embedding dimension: {hidden_size}")

    # ── Load training data ──
    sequences, labels, group_names = load_training_sequences()
    if sequences is None:
        sys.exit(1)

    if args.max_sequences:
        sequences = sequences[:args.max_sequences]
        labels = labels[:args.max_sequences]
        print(f"Limited to {len(sequences)} sequences (--max-sequences)")

    # ── Extract embeddings ──
    embeddings_path = os.path.join(DATA_DIR, "nt_embeddings.npz")

    # Check for cached embeddings
    if os.path.exists(embeddings_path):
        print(f"\nLoading cached embeddings from {embeddings_path}")
        cached = np.load(embeddings_path, allow_pickle=True)
        if (cached["n_sequences"] == len(sequences) and
                str(cached["model_name"]) == model_name):
            X_emb = cached["embeddings"]
            print(f"Loaded {X_emb.shape[0]} cached embeddings ({X_emb.shape[1]}D)")
        else:
            print("Cache mismatch — re-extracting embeddings")
            X_emb = None
    else:
        X_emb = None

    if X_emb is None:
        print(f"\nExtracting NT embeddings for {len(sequences)} sequences...")
        X_emb = np.zeros((len(sequences), hidden_size), dtype=np.float32)

        for i, seq in enumerate(sequences):
            if (i + 1) % 50 == 0 or i == 0:
                print(f"  [{i+1}/{len(sequences)}] "
                      f"({len(seq):,} bp, "
                      f"{max(1, (len(seq) - CHUNK_SIZE) // STRIDE + 1)} chunks)")
            X_emb[i] = extract_nt_embedding(seq, tokenizer, model, device)

        # Cache embeddings
        os.makedirs(DATA_DIR, exist_ok=True)
        np.savez_compressed(
            embeddings_path,
            embeddings=X_emb,
            n_sequences=len(sequences),
            model_name=model_name,
            group_names=np.array(group_names),
        )
        print(f"Cached embeddings: {embeddings_path} "
              f"({os.path.getsize(embeddings_path) / 1024 / 1024:.1f} MB)")

    # ── Train Inc group probe ──
    print(f"\n{'='*60}")
    print("Training Inc group classifier (LogisticRegression probe)")
    print(f"{'='*60}")

    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import cross_val_score
    from sklearn.pipeline import Pipeline

    inc_pipeline = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=2000, C=1.0, multi_class="multinomial",
                                    solver="lbfgs", random_state=42)),
    ])

    # Cross-validation
    scores = cross_val_score(inc_pipeline, X_emb, labels, cv=5, scoring="accuracy")
    print(f"\n5-fold CV accuracy: {scores.mean():.4f} +/- {scores.std():.4f}")

    for i, name in enumerate(group_names):
        n = (labels == i).sum()
        print(f"  {name}: {n} samples")

    # Train on full data
    inc_pipeline.fit(X_emb, labels)

    # Save probe
    inc_probe_path = os.path.join(DATA_DIR, "nt_inc_probe.pkl")
    joblib.dump({
        "pipeline": inc_pipeline,
        "group_names": group_names,
        "model_name": model_name,
        "hidden_size": hidden_size,
        "n_training": len(sequences),
        "cv_accuracy": float(scores.mean()),
    }, inc_probe_path)
    print(f"\nSaved Inc probe: {inc_probe_path}")

    # ── Train AMR probe (if data available) ──
    print(f"\n{'='*60}")
    print("Training AMR class predictor (multi-label probe)")
    print(f"{'='*60}")

    filenames = []
    for inc_name in group_names:
        inc_dir = os.path.join(TRAINING_DIR, inc_name, "fastas")
        for fpath in sorted(glob.glob(os.path.join(inc_dir, "*.fasta"))):
            filenames.append(os.path.basename(fpath))
    if args.max_sequences:
        filenames = filenames[:args.max_sequences]

    amr_labels, amr_classes = load_amr_labels(filenames)

    if amr_labels is not None:
        from sklearn.multioutput import MultiOutputClassifier

        # Filter to classes with sufficient positive examples
        class_counts = amr_labels.sum(axis=0)
        valid_mask = class_counts >= 5
        valid_classes = [c for c, v in zip(amr_classes, valid_mask) if v]
        amr_labels_filtered = amr_labels[:, valid_mask]

        print(f"Training on {len(valid_classes)} AMR classes (>= 5 positive examples)")

        amr_pipeline = Pipeline([
            ("scaler", StandardScaler()),
            ("clf", MultiOutputClassifier(
                LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs",
                                   random_state=42),
                n_jobs=-1,
            )),
        ])

        amr_pipeline.fit(X_emb, amr_labels_filtered)

        amr_probe_path = os.path.join(DATA_DIR, "nt_amr_probe.pkl")
        joblib.dump({
            "pipeline": amr_pipeline,
            "amr_classes": valid_classes,
            "model_name": model_name,
            "hidden_size": hidden_size,
            "n_training": len(sequences),
        }, amr_probe_path)
        print(f"Saved AMR probe: {amr_probe_path}")
        print(f"AMR classes: {valid_classes}")
    else:
        print("Skipping AMR probe (no labeled data available)")

    # ── Summary ──
    print(f"\n{'='*60}")
    print("Training complete!")
    print(f"{'='*60}")
    print(f"Model: {model_name}")
    print(f"Embedding dim: {hidden_size}")
    print(f"Training sequences: {len(sequences)}")
    print(f"Inc group CV accuracy: {scores.mean():.4f}")
    print(f"\nFiles saved:")
    print(f"  {inc_probe_path}")
    if amr_labels is not None:
        print(f"  {os.path.join(DATA_DIR, 'nt_amr_probe.pkl')}")
    print(f"  {embeddings_path}")
    print(f"\nTo use in pLIN app: enable 'Use Nucleotide Transformer' checkbox")


if __name__ == "__main__":
    main()
