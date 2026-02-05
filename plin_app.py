#!/usr/bin/env python3
"""
pLIN Classifier — Streamlit GUI Application
Plasmid Life Identification Number system with AMRFinderPlus integration.
Run: streamlit run plin_app.py

Copyright (C) 2025 Basil Xavier Britto
Licensed under GPL-3.0 with mandatory citation clause.
See LICENSE and CITATION.cff for details.

CITATION REQUIRED: Any use of this software in publications or derivative
works must cite:
    Xavier, B. (2025). pLIN: A Plasmid Life Identification Number System
    for Hierarchical, Permanent Classification of Bacterial Plasmids
    Integrated with Antimicrobial Resistance Gene Surveillance.
    https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification
"""

import os
import io
import glob
import hashlib
import tempfile
import subprocess
import zipfile
import json
import requests
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap
import plotly.graph_objects as go
import plotly.express as px
from itertools import product as iter_product
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from Bio import SeqIO
from sklearn.neighbors import KNeighborsClassifier
from collections import Counter

# ── Page Configuration ────────────────────────────────────────────────────────

st.set_page_config(
    page_title="pLIN Classifier",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Constants ─────────────────────────────────────────────────────────────────

PLIN_THRESHOLDS = {
    "A": 0.150, "B": 0.100, "C": 0.050,
    "D": 0.020, "E": 0.010, "F": 0.001,
}

PLIN_LEVEL_NAMES = {
    "A": "Family", "B": "Subfamily", "C": "Cluster",
    "D": "Subcluster", "E": "Clone", "F": "Strain",
}

ANI_EQUIV = {
    "A": "~85%", "B": "~90%", "C": "~95%",
    "D": "~98%", "E": "~99%", "F": "~99.9%",
}

THRESHOLD_COLORS = {
    "A": "#E53935", "B": "#FB8C00", "C": "#FDD835",
    "D": "#43A047", "E": "#1E88E5", "F": "#8E24AA",
}

STRAIN_COLORS = [
    "#2196F3", "#FF9800", "#4CAF50", "#E91E63", "#9C27B0",
    "#00BCD4", "#FF5722", "#795548", "#607D8B", "#CDDC39",
    "#F44336", "#3F51B5", "#009688", "#FFC107", "#8BC34A",
]

TYPE_COLORS = {"AMR": "#E53935", "STRESS": "#FB8C00", "VIRULENCE": "#8E24AA"}

INC_GROUPS = ["Auto-detect", "IncFII", "IncN", "IncX1", "IncX", "IncH", "Other"]

LINKAGE_METHODS = ["single", "complete", "average", "weighted"]

# Confidence threshold for Inc group classification
# Below this threshold, plasmids are flagged as "Unknown/Novel" Inc type
INC_CONFIDENCE_THRESHOLD = 0.40  # 40% confidence minimum

# Mobility/conjugation marker genes detectable from AMRFinderPlus output
MOBILITY_GENES = {
    "conjugative": {
        "tra": "Transfer (conjugation)",
        "trb": "Transfer (type IV secretion)",
        "vir": "Virulence/T4SS (conjugation-related)",
    },
    "mobilizable": {
        "mob": "Mobilization protein",
        "oriT": "Origin of transfer",
        "nic": "Nickase (relaxase)",
    },
}

# All mobility gene prefixes for scanning
MOBILITY_PREFIXES_CONJUGATIVE = ["traA", "traB", "traC", "traD", "traE", "traF",
                                  "traG", "traH", "traI", "traJ", "traK", "traL",
                                  "traM", "traN", "traP", "traQ", "traR", "traS",
                                  "traT", "traU", "traV", "traW", "traX", "traY",
                                  "trbA", "trbB", "trbC", "trbD", "trbE", "trbF",
                                  "trbG", "trbH", "trbI", "trbJ"]
MOBILITY_PREFIXES_MOBILIZABLE = ["mobA", "mobB", "mobC", "mobD", "mobE", "mobF",
                                  "nikA", "nikB"]

# Paths to precomputed Inc-group classifier data
_APP_DIR = os.path.dirname(os.path.abspath(__file__))
CLASSIFIER_PATH = os.path.join(_APP_DIR, "data", "inc_classifier.npz")
CENTROID_PATH = os.path.join(_APP_DIR, "data", "inc_centroids.npz")

# ── Nucleotide Transformer (optional LLM) ────────────────────────────────────

NT_AVAILABLE = False
try:
    import torch as _torch
    from transformers import AutoTokenizer as _AT, AutoModelForMaskedLM as _AM
    NT_AVAILABLE = True
except ImportError:
    pass

NT_MODELS = {
    "NT-v2-50M (Fast)": "InstaDeepAI/nucleotide-transformer-v2-50m-multi-species",
    "NT-v2-100M": "InstaDeepAI/nucleotide-transformer-v2-100m-multi-species",
    "NT-v2-250M": "InstaDeepAI/nucleotide-transformer-v2-250m-multi-species",
    "NT-v2-500M (Best)": "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
}
NT_INC_PROBE_PATH = os.path.join(_APP_DIR, "data", "nt_inc_probe.pkl")
NT_AMR_PROBE_PATH = os.path.join(_APP_DIR, "data", "nt_amr_probe.pkl")
NT_CHUNK_SIZE = 5000
NT_STRIDE = 2500


# ══════════════════════════════════════════════════════════════════════════════
#  INC GROUP CLASSIFICATION
# ══════════════════════════════════════════════════════════════════════════════

@st.cache_resource(show_spinner=False)
def load_inc_classifier():
    """Load precomputed Inc-group KNN classifier (96% CV accuracy)."""
    if os.path.exists(CLASSIFIER_PATH):
        data = np.load(CLASSIFIER_PATH, allow_pickle=True)
        X_train = data["X"]
        y_train = data["y"]
        group_names = [str(g) for g in data["group_names"]]
        knn = KNeighborsClassifier(n_neighbors=5, metric="cosine", weights="distance")
        knn.fit(X_train, y_train)
        return group_names, knn
    elif os.path.exists(CENTROID_PATH):
        data = np.load(CENTROID_PATH, allow_pickle=True)
        group_names = [str(g) for g in data["group_names"]]
        return group_names, data["centroids"]
    return None, None


def _kmer_vector_single(sequence):
    """Compute normalised 4-mer frequency vector for a single sequence."""
    bases = "ACGT"
    all_kmers = ["".join(p) for p in iter_product(bases, repeat=4)]
    kmer_idx = {km: i for i, km in enumerate(all_kmers)}
    seq = sequence.upper()
    counts = np.zeros(256, dtype=np.float64)
    for i in range(len(seq) - 3):
        kmer = seq[i:i + 4]
        if kmer in kmer_idx:
            counts[kmer_idx[kmer]] += 1
    total = counts.sum()
    if total > 0:
        counts /= total
    return counts


def classify_inc_group(sequence, group_names, classifier):
    """Classify a plasmid to its Inc group using KNN or centroid distance.

    Returns dict with:
        - predicted_group: str (best Inc type or "Unknown/Novel" if low confidence)
        - confidence: float (0-1)
        - is_low_confidence: bool
        - top5_candidates: list of (inc_type, confidence) tuples
        - all_probabilities: dict {inc_type: probability}
    """
    vec = _kmer_vector_single(sequence).reshape(1, -1).astype(np.float32)

    if isinstance(classifier, KNeighborsClassifier):
        proba = classifier.predict_proba(vec)[0]
        pred_idx = np.argmax(proba)
        best_group = group_names[pred_idx]
        confidence = float(proba[pred_idx])
        proba_dict = {group_names[j]: round(float(proba[j]), 4) for j in range(len(group_names))}

        # Get top 5 candidates sorted by confidence
        sorted_candidates = sorted(proba_dict.items(), key=lambda x: x[1], reverse=True)[:5]

        # Check if confidence is below threshold
        is_low_confidence = confidence < INC_CONFIDENCE_THRESHOLD

        return {
            "predicted_group": "Unknown/Novel" if is_low_confidence else best_group,
            "best_match": best_group,  # Always store the best match even if flagged as Unknown
            "confidence": confidence,
            "is_low_confidence": is_low_confidence,
            "top5_candidates": sorted_candidates,
            "all_probabilities": proba_dict,
        }
    else:
        # Centroid fallback — use cosine similarity as confidence
        centroids = classifier
        v = vec.flatten()
        similarities = {}
        for j, name in enumerate(group_names):
            c = centroids[j]
            norm_a, norm_b = np.linalg.norm(v), np.linalg.norm(c)
            sim = np.dot(v, c) / (norm_a * norm_b) if norm_a > 0 and norm_b > 0 else 0.0
            similarities[name] = round(sim, 4)
        best_group = max(similarities, key=similarities.get)
        confidence = similarities[best_group]

        # Get top 5 candidates sorted by similarity
        sorted_candidates = sorted(similarities.items(), key=lambda x: x[1], reverse=True)[:5]

        # Check if confidence is below threshold
        is_low_confidence = confidence < INC_CONFIDENCE_THRESHOLD

        return {
            "predicted_group": "Unknown/Novel" if is_low_confidence else best_group,
            "best_match": best_group,
            "confidence": confidence,
            "is_low_confidence": is_low_confidence,
            "top5_candidates": sorted_candidates,
            "all_probabilities": similarities,
        }


# ══════════════════════════════════════════════════════════════════════════════
#  ADAPTIVE THRESHOLD CALIBRATION
# ══════════════════════════════════════════════════════════════════════════════

@st.cache_data(show_spinner=False)
def calibrate_inc_thresholds():
    """Calibrate pLIN thresholds per Inc group from training data distance distributions.

    Uses quantile-based calibration: for each Inc group, compute pairwise cosine
    distances and set thresholds at specific quantile percentiles that correspond
    to the hierarchical levels (Family→Strain).

    Returns dict: {inc_group: {A: thresh, B: thresh, ...}} or None if no data.
    """
    if not os.path.exists(CLASSIFIER_PATH):
        return None

    data = np.load(CLASSIFIER_PATH, allow_pickle=True)
    X = data["X"]
    y = data["y"]
    group_names = [str(g) for g in data["group_names"]]

    # Quantile percentiles for each pLIN level (from broadest to finest)
    # These correspond to the fraction of within-group distances that should
    # fall below each threshold
    level_quantiles = {
        "A": 0.99,   # Family — nearly all within-group distances below this
        "B": 0.95,   # Subfamily
        "C": 0.75,   # Cluster
        "D": 0.50,   # Subcluster — median distance
        "E": 0.25,   # Clone
        "F": 0.05,   # Strain — only very close pairs
    }

    calibrated = {}
    for i, name in enumerate(group_names):
        mask = y == i
        X_group = X[mask]
        if X_group.shape[0] < 10:
            continue

        # Sample to keep computation tractable
        n = X_group.shape[0]
        if n > 500:
            rng = np.random.default_rng(42)
            idx = rng.choice(n, 500, replace=False)
            X_sub = X_group[idx]
        else:
            X_sub = X_group

        dists = pdist(X_sub.astype(np.float64), metric="cosine")

        thresholds = {}
        for level, q in level_quantiles.items():
            val = float(np.quantile(dists, q))
            # Ensure minimum separation between levels
            thresholds[level] = round(max(val, 0.0005), 6)

        # Ensure monotonic: A > B > C > D > E > F
        levels = list("ABCDEF")
        for j in range(1, len(levels)):
            if thresholds[levels[j]] >= thresholds[levels[j-1]]:
                thresholds[levels[j]] = round(thresholds[levels[j-1]] * 0.7, 6)

        calibrated[name] = thresholds

    return calibrated if calibrated else None


# ══════════════════════════════════════════════════════════════════════════════
#  MOBILITY PREDICTION
# ══════════════════════════════════════════════════════════════════════════════

def classify_mobility(amr_df, source_file):
    """Classify plasmid mobility from AMRFinderPlus output.

    Categories:
    - Conjugative: has tra/trb transfer genes (can self-transfer)
    - Mobilizable: has mob genes but no full tra system (needs helper)
    - Non-mobilizable: no detectable transfer/mobilization genes

    Returns (category, genes_found, details).
    """
    if amr_df is None or len(amr_df) == 0:
        return "Unknown", [], "No AMR data"

    sf = source_file.replace(".fasta", "").replace(".fa", "").replace(".fna", "")
    hits = amr_df[amr_df["source_file"] == sf]

    if len(hits) == 0:
        return "Non-mobilizable", [], "No genes detected"

    gene_col = "Element symbol" if "Element symbol" in hits.columns else None
    name_col = "Element name" if "Element name" in hits.columns else None

    if gene_col is None:
        return "Unknown", [], "No gene symbol column"

    all_genes = hits[gene_col].tolist()
    all_names = hits[name_col].tolist() if name_col else []

    # Check for conjugative markers
    conj_found = []
    for gene in all_genes:
        for prefix in MOBILITY_PREFIXES_CONJUGATIVE:
            if gene.startswith(prefix):
                conj_found.append(gene)
                break

    # Also check gene names/descriptions for transfer keywords
    for name in all_names:
        name_lower = str(name).lower()
        if "conjugal transfer" in name_lower or "type iv secretion" in name_lower:
            # Find the corresponding gene
            idx = all_names.index(name)
            g = all_genes[idx]
            if g not in conj_found:
                conj_found.append(g)

    # Check for mobilizable markers
    mob_found = []
    for gene in all_genes:
        for prefix in MOBILITY_PREFIXES_MOBILIZABLE:
            if gene.startswith(prefix):
                mob_found.append(gene)
                break

    if conj_found:
        return "Conjugative", list(set(conj_found)), f"{len(set(conj_found))} transfer gene(s)"
    elif mob_found:
        return "Mobilizable", list(set(mob_found)), f"{len(set(mob_found))} mobilization gene(s)"
    else:
        return "Non-mobilizable", [], "No transfer/mobilization genes detected"


# ══════════════════════════════════════════════════════════════════════════════
#  OUTBREAK / CLONE DETECTION
# ══════════════════════════════════════════════════════════════════════════════

def detect_outbreak_clusters(plin_df, integrated_df):
    """Flag potential outbreak clusters: plasmids sharing same pLIN strain code (F)
    AND identical AMR resistance profile.

    Returns list of dicts with cluster info.
    """
    if integrated_df is None or len(integrated_df) == 0:
        return []

    df = integrated_df.copy()

    # Create AMR fingerprint: sorted set of AMR genes
    if "AMR_genes" in df.columns:
        df["_amr_fingerprint"] = df["AMR_genes"].fillna("").apply(
            lambda x: "|".join(sorted(g.strip() for g in x.split(";") if g.strip()))
        )
    else:
        df["_amr_fingerprint"] = ""

    # Group by strain cluster (bin_F) + AMR fingerprint
    if "bin_F" not in df.columns:
        # Extract from pLIN code
        df["bin_F"] = df["pLIN"].apply(lambda x: x.split(".")[-1] if isinstance(x, str) else "")

    clusters = []
    grouped = df.groupby(["bin_F", "_amr_fingerprint"])
    for (strain_f, amr_fp), group in grouped:
        if len(group) >= 2 and amr_fp:  # At least 2 plasmids with same AMR
            plasmids = group["plasmid_id"].tolist()
            amr_genes = [g.strip() for g in amr_fp.split("|") if g.strip()]
            plin_code = group["pLIN"].iloc[0]
            clusters.append({
                "strain_cluster": int(strain_f) if str(strain_f).isdigit() else strain_f,
                "pLIN": plin_code,
                "n_plasmids": len(group),
                "plasmids": plasmids,
                "amr_genes": amr_genes,
                "n_amr_genes": len(amr_genes),
                "risk_level": "HIGH" if len(amr_genes) >= 3 else "MODERATE",
            })

    # Sort by risk
    clusters.sort(key=lambda c: (-c["n_amr_genes"], -c["n_plasmids"]))
    return clusters


# ══════════════════════════════════════════════════════════════════════════════
#  NUCLEOTIDE TRANSFORMER (OPTIONAL LLM)
# ══════════════════════════════════════════════════════════════════════════════

@st.cache_resource(show_spinner="Loading Nucleotide Transformer model...")
def load_nt_model(model_name, device_name="cpu"):
    """Load pre-trained Nucleotide Transformer model and tokenizer."""
    import torch
    from transformers import AutoTokenizer, AutoModelForMaskedLM

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForMaskedLM.from_pretrained(model_name)
    model.eval()
    device = torch.device(device_name)
    model.to(device)
    return tokenizer, model, device


def extract_nt_embedding(sequence, tokenizer, model, device,
                         chunk_size=NT_CHUNK_SIZE, stride=NT_STRIDE):
    """Extract mean-pooled NT embedding for a single plasmid sequence.

    Long sequences are split into overlapping chunks, embedded independently,
    and chunk embeddings are averaged to produce a single vector.
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
        hidden = outputs.hidden_states[-1]
        mask = tokens["attention_mask"].unsqueeze(-1).float()
        mean_emb = (hidden * mask).sum(dim=1) / mask.sum(dim=1)
        embeddings.append(mean_emb.squeeze().cpu().numpy())

    return np.mean(embeddings, axis=0).astype(np.float32)


def run_nt_predictions(sequences, tokenizer, model, device, progress_cb=None):
    """Extract NT embeddings and run probe predictions for all sequences.

    Returns dict with embeddings, inc predictions, and amr predictions.
    """
    import joblib

    # Extract embeddings
    embeddings = []
    for i, seq in enumerate(sequences):
        emb = extract_nt_embedding(seq, tokenizer, model, device)
        embeddings.append(emb)
        if progress_cb:
            progress_cb((i + 1) / len(sequences))
    X = np.array(embeddings, dtype=np.float32)

    results = {"embeddings": X}

    # Inc group probe
    if os.path.exists(NT_INC_PROBE_PATH):
        probe_data = joblib.load(NT_INC_PROBE_PATH)
        pipeline = probe_data["pipeline"]
        group_names = probe_data["group_names"]
        inc_preds = pipeline.predict(X)
        inc_proba = pipeline.predict_proba(X)
        results["inc_preds"] = [group_names[p] for p in inc_preds]
        results["inc_proba"] = inc_proba
        results["inc_groups"] = group_names
        results["inc_cv_accuracy"] = probe_data.get("cv_accuracy", None)
    else:
        results["inc_preds"] = None

    # AMR class probe
    if os.path.exists(NT_AMR_PROBE_PATH):
        amr_data = joblib.load(NT_AMR_PROBE_PATH)
        amr_pipeline = amr_data["pipeline"]
        amr_classes = amr_data["amr_classes"]
        amr_preds = amr_pipeline.predict(X)
        results["amr_preds"] = amr_preds
        results["amr_classes"] = amr_classes
    else:
        results["amr_preds"] = None

    return results


def detect_nt_device():
    """Auto-detect best available compute device for NT inference."""
    if not NT_AVAILABLE:
        return "cpu"
    import torch
    if torch.cuda.is_available():
        return "cuda"
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return "mps"
    return "cpu"


# ══════════════════════════════════════════════════════════════════════════════
#  CORE PIPELINE FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def parse_uploaded_fastas(uploaded_files, inc_type):
    """Parse uploaded FASTA files. Auto-detects Inc group when inc_type='Auto-detect'."""
    auto_detect = (inc_type == "Auto-detect")
    group_names, classifier = None, None
    if auto_detect:
        group_names, classifier = load_inc_classifier()
        if group_names is None:
            auto_detect = False
            inc_type = "Unknown"

    records = []
    for uf in uploaded_files:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", mode="wb") as tmp:
            tmp.write(uf.getvalue())
            tmp_path = tmp.name
        try:
            for rec in SeqIO.parse(tmp_path, "fasta"):
                seq = str(rec.seq)
                rec_dict = {
                    "plasmid_id": rec.id,
                    "sequence": seq,
                    "length": len(seq),
                    "source_file": uf.name,
                }
                if auto_detect:
                    result = classify_inc_group(seq, group_names, classifier)
                    rec_dict["inc_type"] = result["predicted_group"]
                    rec_dict["inc_best_match"] = result["best_match"]
                    rec_dict["inc_confidence"] = round(result["confidence"], 4)
                    rec_dict["inc_is_low_confidence"] = result["is_low_confidence"]
                    rec_dict["inc_top5_candidates"] = result["top5_candidates"]
                    rec_dict["inc_probabilities"] = result["all_probabilities"]
                else:
                    rec_dict["inc_type"] = inc_type

                records.append(rec_dict)
        finally:
            os.unlink(tmp_path)
    return records


@st.cache_data(show_spinner=False)
def compute_kmer_vectors(_sequences, k=4):
    """Compute normalised tetranucleotide frequency vectors."""
    bases = "ACGT"
    all_kmers = ["".join(p) for p in iter_product(bases, repeat=k)]
    vectors = np.zeros((len(_sequences), len(all_kmers)), dtype=np.float64)
    for idx, seq in enumerate(_sequences):
        s = seq.upper()
        total = max(len(s) - k + 1, 1)
        for ki, kmer in enumerate(all_kmers):
            vectors[idx, ki] = s.count(kmer) / total
    return vectors


def assign_plin_codes(vectors, linkage_method="single", thresholds=None):
    """Cluster plasmids and assign hierarchical pLIN codes.

    Args:
        vectors: k-mer frequency matrix
        linkage_method: 'single', 'complete', 'average', or 'weighted'
        thresholds: dict {A: val, ...} or None to use defaults
    """
    active_thresholds = thresholds if thresholds else PLIN_THRESHOLDS
    dist_condensed = pdist(vectors, metric="cosine")
    Z = linkage(dist_condensed, method=linkage_method)

    cluster_assignments = {}
    for bname, thresh in active_thresholds.items():
        cluster_assignments[bname] = fcluster(Z, t=thresh, criterion="distance")

    n = vectors.shape[0]
    plin_codes = []
    for i in range(n):
        parts = [str(cluster_assignments[b][i]) for b in active_thresholds]
        plin_codes.append(".".join(parts))

    return plin_codes, cluster_assignments, Z, dist_condensed


def build_results_df(records, plin_codes, cluster_assignments):
    """Build results DataFrame."""
    rows = []
    for i, rec in enumerate(records):
        row = {
            "plasmid_id": rec["plasmid_id"],
            "inc_type": rec["inc_type"],
            "source_file": rec["source_file"],
            "length_bp": rec["length"],
            "pLIN": plin_codes[i],
        }
        if "inc_confidence" in rec:
            row["inc_confidence"] = rec["inc_confidence"]
        if "inc_is_low_confidence" in rec:
            row["inc_is_low_confidence"] = rec["inc_is_low_confidence"]
        if "inc_best_match" in rec:
            row["inc_best_match"] = rec["inc_best_match"]
        if "inc_top5_candidates" in rec:
            # Format top 5 as string for display: "IncFII (45%), IncN (30%), ..."
            top5_str = ", ".join([f"{inc} ({conf*100:.1f}%)" for inc, conf in rec["inc_top5_candidates"]])
            row["inc_top5_candidates"] = top5_str
        for b in PLIN_THRESHOLDS:
            row[f"bin_{b}"] = int(cluster_assignments[b][i])
        rows.append(row)
    return pd.DataFrame(rows)


# ══════════════════════════════════════════════════════════════════════════════
#  AMRFINDERPLUS FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def detect_amrfinder():
    """Auto-detect AMRFinderPlus binary and database."""
    binary = None
    database = None

    # 1. Check PATH
    try:
        result = subprocess.run(["which", "amrfinder"], capture_output=True, text=True)
        if result.returncode == 0 and result.stdout.strip():
            binary = result.stdout.strip()
    except Exception:
        pass

    # 2. Check conda envs
    if not binary:
        home = os.path.expanduser("~")
        search_dirs = [
            os.path.join(home, "miniconda3", "envs"),
            os.path.join(home, "miniforge3", "envs"),
            os.path.join(home, "anaconda3", "envs"),
            os.path.join(home, "mambaforge", "envs"),
        ]
        for base in search_dirs:
            if os.path.isdir(base):
                for env in sorted(os.listdir(base)):
                    candidate = os.path.join(base, env, "bin", "amrfinder")
                    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                        binary = candidate
                        break
            if binary:
                break

    # Find database
    if binary:
        prefix = os.path.dirname(os.path.dirname(binary))
        db_base = os.path.join(prefix, "share", "amrfinderplus", "data")
        if os.path.isdir(db_base):
            versions = sorted(
                [d for d in os.listdir(db_base) if d.startswith("20") and os.path.isdir(os.path.join(db_base, d))]
            )
            if versions:
                database = os.path.join(db_base, versions[-1])

    return binary, database


def run_amrfinder_on_files(uploaded_files, binary, database, progress_callback=None):
    """Run AMRFinderPlus on uploaded FASTA files."""
    all_results = []
    total = len(uploaded_files)

    with tempfile.TemporaryDirectory() as tmpdir:
        for idx, uf in enumerate(uploaded_files):
            fasta_path = os.path.join(tmpdir, uf.name)
            with open(fasta_path, "wb") as f:
                f.write(uf.getvalue())

            out_path = os.path.join(tmpdir, f"{uf.name}.amr.tsv")

            cmd = [binary, "-n", fasta_path, "--plus", "-o", out_path]
            if database:
                cmd.extend(["-d", database])

            try:
                subprocess.run(cmd, capture_output=True, timeout=300)
                if os.path.isfile(out_path):
                    df = pd.read_csv(out_path, sep="\t")
                    df.insert(0, "source_file", uf.name.replace(".fasta", "").replace(".fa", "").replace(".fna", ""))
                    all_results.append(df)
            except Exception:
                pass

            if progress_callback:
                progress_callback((idx + 1) / total)

    if all_results:
        return pd.concat(all_results, ignore_index=True)
    return pd.DataFrame()


def integrate_plin_amr(plin_df, amr_df):
    """Integrate pLIN assignments with AMR results."""
    rows = []
    for _, row in plin_df.iterrows():
        sf = row["source_file"].replace(".fasta", "").replace(".fa", "").replace(".fna", "")
        hits = amr_df[amr_df["source_file"] == sf] if len(amr_df) > 0 else pd.DataFrame()

        amr_hits = hits[hits["Type"] == "AMR"] if "Type" in hits.columns and len(hits) > 0 else pd.DataFrame()
        stress_hits = hits[hits["Type"] == "STRESS"] if "Type" in hits.columns and len(hits) > 0 else pd.DataFrame()
        vir_hits = hits[hits["Type"] == "VIRULENCE"] if "Type" in hits.columns and len(hits) > 0 else pd.DataFrame()

        rows.append({
            **row.to_dict(),
            "total_hits": len(hits),
            "AMR_count": len(amr_hits),
            "STRESS_count": len(stress_hits),
            "VIR_count": len(vir_hits),
            "AMR_genes": "; ".join(sorted(amr_hits["Element symbol"].unique())) if len(amr_hits) > 0 else "",
            "STRESS_genes": "; ".join(sorted(stress_hits["Element symbol"].unique())) if len(stress_hits) > 0 else "",
            "AMR_classes": "; ".join(sorted(amr_hits["Class"].unique())) if len(amr_hits) > 0 else "",
        })
    return pd.DataFrame(rows)


# ══════════════════════════════════════════════════════════════════════════════
#  PRODIGAL GENE ANNOTATION FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def detect_prodigal():
    """Auto-detect Prodigal binary."""
    binary = None

    # 1. Check PATH
    try:
        result = subprocess.run(["which", "prodigal"], capture_output=True, text=True)
        if result.returncode == 0 and result.stdout.strip():
            binary = result.stdout.strip()
    except Exception:
        pass

    # 2. Check conda envs
    if not binary:
        home = os.path.expanduser("~")
        search_dirs = [
            os.path.join(home, "miniconda3", "envs"),
            os.path.join(home, "miniforge3", "envs"),
            os.path.join(home, "anaconda3", "envs"),
            os.path.join(home, "mambaforge", "envs"),
        ]
        for base in search_dirs:
            if os.path.isdir(base):
                for env in sorted(os.listdir(base)):
                    candidate = os.path.join(base, env, "bin", "prodigal")
                    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                        binary = candidate
                        break
            if binary:
                break

    return binary


def parse_prodigal_gff(gff_path):
    """Parse Prodigal GFF output to extract gene information."""
    genes = []
    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "CDS":
                continue

            seqid = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(";"):
                if "=" in attr:
                    key, val = attr.split("=", 1)
                    attr_dict[key] = val

            gene_id = attr_dict.get("ID", f"{seqid}_{start}_{end}")
            partial = attr_dict.get("partial", "00")

            genes.append({
                "gene_id": gene_id,
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strand,
                "length_bp": end - start + 1,
                "length_aa": (end - start + 1) // 3,
                "partial": partial,
                "is_complete": partial == "00",
            })

    return genes


def run_prodigal_on_files(uploaded_files, binary, progress_callback=None):
    """Run Prodigal on uploaded FASTA files.

    Returns:
        - genes_df: DataFrame with all predicted genes
        - summary_df: DataFrame with per-plasmid summary stats
    """
    all_genes = []
    summaries = []
    total = len(uploaded_files)

    with tempfile.TemporaryDirectory() as tmpdir:
        for idx, uf in enumerate(uploaded_files):
            fasta_path = os.path.join(tmpdir, uf.name)
            with open(fasta_path, "wb") as f:
                f.write(uf.getvalue())

            gff_path = os.path.join(tmpdir, f"{uf.name}.gff")
            proteins_path = os.path.join(tmpdir, f"{uf.name}.faa")

            # Run Prodigal in metagenomic mode (better for plasmids with variable GC)
            cmd = [
                binary,
                "-i", fasta_path,
                "-o", gff_path,
                "-a", proteins_path,
                "-f", "gff",
                "-p", "meta",  # Metagenomic mode - good for diverse plasmids
                "-q",  # Quiet mode
            ]

            source_name = uf.name.replace(".fasta", "").replace(".fa", "").replace(".fna", "")

            try:
                subprocess.run(cmd, capture_output=True, timeout=300)
                if os.path.isfile(gff_path):
                    genes = parse_prodigal_gff(gff_path)

                    # Get sequence length from FASTA
                    seq_length = 0
                    for rec in SeqIO.parse(fasta_path, "fasta"):
                        seq_length += len(rec.seq)

                    # Add source file to each gene
                    for g in genes:
                        g["source_file"] = source_name
                    all_genes.extend(genes)

                    # Calculate summary stats
                    total_genes = len(genes)
                    complete_genes = sum(1 for g in genes if g["is_complete"])
                    total_coding_bp = sum(g["length_bp"] for g in genes)
                    coding_density = (total_coding_bp / seq_length * 100) if seq_length > 0 else 0
                    avg_gene_length = np.mean([g["length_aa"] for g in genes]) if genes else 0

                    summaries.append({
                        "source_file": source_name,
                        "sequence_length": seq_length,
                        "total_genes": total_genes,
                        "complete_genes": complete_genes,
                        "partial_genes": total_genes - complete_genes,
                        "total_coding_bp": total_coding_bp,
                        "coding_density_pct": round(coding_density, 1),
                        "avg_gene_length_aa": round(avg_gene_length, 1),
                    })
                else:
                    # No output - add empty summary
                    summaries.append({
                        "source_file": source_name,
                        "sequence_length": 0,
                        "total_genes": 0,
                        "complete_genes": 0,
                        "partial_genes": 0,
                        "total_coding_bp": 0,
                        "coding_density_pct": 0,
                        "avg_gene_length_aa": 0,
                    })
            except Exception:
                summaries.append({
                    "source_file": source_name,
                    "sequence_length": 0,
                    "total_genes": 0,
                    "complete_genes": 0,
                    "partial_genes": 0,
                    "total_coding_bp": 0,
                    "coding_density_pct": 0,
                    "avg_gene_length_aa": 0,
                })

            if progress_callback:
                progress_callback((idx + 1) / total)

    genes_df = pd.DataFrame(all_genes) if all_genes else pd.DataFrame()
    summary_df = pd.DataFrame(summaries) if summaries else pd.DataFrame()

    return genes_df, summary_df


# ══════════════════════════════════════════════════════════════════════════════
#  BACTERIAL BUDDY — OLLAMA LLM CHATBOT
# ══════════════════════════════════════════════════════════════════════════════

OLLAMA_DEFAULT_URL = "http://localhost:11434"
OLLAMA_MODELS = ["llama3.2", "llama3.1", "llama3", "mistral", "mixtral", "gemma2", "phi3"]

BACTERIAL_BUDDY_SYSTEM_PROMPT = """You are Bacterial Buddy, a friendly and knowledgeable AI assistant specialized in plasmid biology and antimicrobial resistance. You help researchers understand their plasmid analysis results from the pLIN (Plasmid Life Identification Number) classification tool.

Your personality:
- Friendly, approachable, and enthusiastic about microbiology
- Use occasional bacteria-themed humor when appropriate
- Explain complex concepts in accessible terms
- Always be scientifically accurate

Your expertise includes:
- Plasmid biology, replication, and evolution
- Incompatibility (Inc) groups and their significance
- Antimicrobial resistance genes and mechanisms
- Plasmid mobility and horizontal gene transfer
- Epidemiological implications of plasmid spread
- Interpreting pLIN codes and clustering results

When answering questions:
1. Use the analysis context provided to give specific, relevant answers
2. If asked about specific plasmids or genes, refer to the data
3. Explain the clinical/epidemiological significance when relevant
4. Suggest next steps or further analyses when appropriate

Keep responses concise but informative. Use bullet points for lists."""


def detect_ollama():
    """Check if Ollama is running and return available models."""
    try:
        response = requests.get(f"{OLLAMA_DEFAULT_URL}/api/tags", timeout=2)
        if response.status_code == 200:
            data = response.json()
            models = [m["name"].split(":")[0] for m in data.get("models", [])]
            return True, list(set(models))
    except Exception:
        pass
    return False, []


def build_analysis_context(session_state):
    """Build a context string from the current analysis results."""
    context_parts = []

    # Basic analysis info
    plin_df = session_state.get("plin_df")
    if plin_df is not None and len(plin_df) > 0:
        context_parts.append(f"**Analysis Summary:**")
        context_parts.append(f"- Total plasmids analyzed: {len(plin_df)}")
        context_parts.append(f"- Unique pLIN codes: {plin_df['pLIN'].nunique()}")

        # Inc groups
        if "inc_type" in plin_df.columns:
            inc_counts = plin_df["inc_type"].value_counts().to_dict()
            inc_str = ", ".join([f"{k}: {v}" for k, v in inc_counts.items()])
            context_parts.append(f"- Inc groups detected: {inc_str}")

        # Low confidence warnings
        if "inc_is_low_confidence" in plin_df.columns:
            low_conf = plin_df["inc_is_low_confidence"].sum()
            if low_conf > 0:
                context_parts.append(f"- Low-confidence Inc predictions (Unknown/Novel): {low_conf}")

    # AMR data
    amr_df = session_state.get("amr_df")
    if amr_df is not None and len(amr_df) > 0:
        context_parts.append(f"\n**AMR Analysis:**")
        context_parts.append(f"- Total AMR/stress/virulence hits: {len(amr_df)}")

        if "Type" in amr_df.columns:
            type_counts = amr_df["Type"].value_counts().to_dict()
            context_parts.append(f"- AMR genes: {type_counts.get('AMR', 0)}")
            context_parts.append(f"- Stress genes: {type_counts.get('STRESS', 0)}")
            context_parts.append(f"- Virulence genes: {type_counts.get('VIRULENCE', 0)}")

        if "Class" in amr_df.columns:
            top_classes = amr_df["Class"].value_counts().head(5).to_dict()
            classes_str = ", ".join([f"{k}: {v}" for k, v in top_classes.items()])
            context_parts.append(f"- Top AMR drug classes: {classes_str}")

        if "Element symbol" in amr_df.columns:
            top_genes = amr_df["Element symbol"].value_counts().head(10).to_dict()
            genes_str = ", ".join([f"{k} ({v})" for k, v in top_genes.items()])
            context_parts.append(f"- Most common genes: {genes_str}")

    # Mobility
    mob_df = session_state.get("mobility_results")
    if mob_df is not None and len(mob_df) > 0:
        context_parts.append(f"\n**Mobility Prediction:**")
        mob_counts = mob_df["mobility"].value_counts().to_dict()
        context_parts.append(f"- Conjugative: {mob_counts.get('Conjugative', 0)}")
        context_parts.append(f"- Mobilizable: {mob_counts.get('Mobilizable', 0)}")
        context_parts.append(f"- Non-mobilizable: {mob_counts.get('Non-mobilizable', 0)}")

    # Outbreak clusters
    outbreak = session_state.get("outbreak_clusters", [])
    if outbreak:
        context_parts.append(f"\n**Outbreak Detection:**")
        context_parts.append(f"- Potential outbreak clusters detected: {len(outbreak)}")
        for i, cluster in enumerate(outbreak[:3]):  # Show top 3
            context_parts.append(f"  - Cluster {i+1}: {cluster.get('count', '?')} plasmids, "
                               f"risk: {cluster.get('risk_level', '?')}")

    # Prodigal
    prodigal_sum = session_state.get("prodigal_summary_df")
    if prodigal_sum is not None and len(prodigal_sum) > 0:
        context_parts.append(f"\n**Gene Annotation (Prodigal):**")
        total_genes = prodigal_sum["total_genes"].sum()
        avg_density = prodigal_sum["coding_density_pct"].mean()
        context_parts.append(f"- Total predicted genes: {total_genes}")
        context_parts.append(f"- Average coding density: {avg_density:.1f}%")

    # Plasmid list
    if plin_df is not None and len(plin_df) > 0 and len(plin_df) <= 20:
        context_parts.append(f"\n**Plasmid Details:**")
        for _, row in plin_df.iterrows():
            detail = f"- {row['plasmid_id']}: pLIN {row['pLIN']}"
            if "inc_type" in row:
                detail += f", Inc: {row['inc_type']}"
            context_parts.append(detail)

    return "\n".join(context_parts) if context_parts else "No analysis has been run yet."


def chat_with_ollama(messages, model, context=""):
    """Send messages to Ollama and get a response."""
    # Build the system message with context
    system_msg = BACTERIAL_BUDDY_SYSTEM_PROMPT
    if context:
        system_msg += f"\n\n**Current Analysis Context:**\n{context}"

    # Prepare messages for Ollama
    ollama_messages = [{"role": "system", "content": system_msg}]
    for msg in messages:
        ollama_messages.append({"role": msg["role"], "content": msg["content"]})

    try:
        response = requests.post(
            f"{OLLAMA_DEFAULT_URL}/api/chat",
            json={
                "model": model,
                "messages": ollama_messages,
                "stream": False,
            },
            timeout=60,
        )
        if response.status_code == 200:
            data = response.json()
            return data.get("message", {}).get("content", "I couldn't generate a response.")
        else:
            return f"Error: Ollama returned status {response.status_code}"
    except requests.exceptions.Timeout:
        return "The response took too long. Try a simpler question or a smaller model."
    except Exception as e:
        return f"Error connecting to Ollama: {str(e)}"


def stream_chat_with_ollama(messages, model, context=""):
    """Stream responses from Ollama for a better UX."""
    system_msg = BACTERIAL_BUDDY_SYSTEM_PROMPT
    if context:
        system_msg += f"\n\n**Current Analysis Context:**\n{context}"

    ollama_messages = [{"role": "system", "content": system_msg}]
    for msg in messages:
        ollama_messages.append({"role": msg["role"], "content": msg["content"]})

    try:
        response = requests.post(
            f"{OLLAMA_DEFAULT_URL}/api/chat",
            json={
                "model": model,
                "messages": ollama_messages,
                "stream": True,
            },
            stream=True,
            timeout=120,
        )
        if response.status_code == 200:
            for line in response.iter_lines():
                if line:
                    data = json.loads(line)
                    if "message" in data and "content" in data["message"]:
                        yield data["message"]["content"]
                    if data.get("done", False):
                        break
        else:
            yield f"Error: Ollama returned status {response.status_code}"
    except Exception as e:
        yield f"Error: {str(e)}"


# ══════════════════════════════════════════════════════════════════════════════
#  VISUALIZATION FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def get_strain_cmap(strain_clusters):
    unique_strains = sorted(set(strain_clusters))
    return {s: STRAIN_COLORS[i % len(STRAIN_COLORS)] for i, s in enumerate(unique_strains)}


def fig_to_bytes(fig, fmt="png", dpi=300):
    """Convert matplotlib figure to bytes buffer."""
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=dpi, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    return buf.getvalue()


def plot_rectangular_cladogram(Z, labels, plin_codes, strain_clusters):
    strain_cmap = get_strain_cmap(strain_clusters)
    unique_strains = sorted(set(strain_clusters))

    fig, (ax_dendro, ax_plin) = plt.subplots(
        1, 2, figsize=(14, max(8, len(labels) * 0.4)),
        gridspec_kw={"width_ratios": [3, 2], "wspace": 0.02},
    )

    with plt.rc_context({"lines.linewidth": 2.0}):
        ddata = dendrogram(Z, labels=labels, orientation="right",
                           leaf_font_size=9, ax=ax_dendro,
                           color_threshold=0, above_threshold_color="#444444")

    max_merge = Z[:, 2].max()
    ax_dendro.set_xlim(0, max_merge * 1.15)

    leaf_colors = {labels[i]: strain_cmap[sc] for i, sc in enumerate(strain_clusters)}
    for lbl in ax_dendro.get_yticklabels():
        txt = lbl.get_text()
        if txt in leaf_colors:
            lbl.set_color(leaf_colors[txt])
            lbl.set_fontweight("bold")

    for bname, thresh in PLIN_THRESHOLDS.items():
        if thresh <= max_merge * 1.15:
            ax_dendro.axvline(x=thresh, color=THRESHOLD_COLORS[bname],
                              linestyle="--", linewidth=1.2, alpha=0.7)

    offscreen = [(b, t) for b, t in PLIN_THRESHOLDS.items() if t > max_merge * 1.15]
    if offscreen:
        lines = ["Thresholds beyond range:"] + [f"  {b} ({PLIN_LEVEL_NAMES[b]}): d\u2264{t}" for b, t in offscreen]
        ax_dendro.text(0.98, 0.02, "\n".join(lines), transform=ax_dendro.transAxes,
                       fontsize=7, va="bottom", ha="right",
                       bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.9))

    ax_dendro.set_xlabel("Cosine Distance (4-mer)", fontsize=10, fontweight="bold")
    ax_dendro.set_title("Dendrogram", fontsize=11, fontweight="bold")

    ax_plin.set_ylim(ax_dendro.get_ylim())
    ax_plin.set_xlim(0, 1)
    ax_plin.axis("off")
    ax_plin.set_title("pLIN Code", fontsize=11, fontweight="bold")

    leaf_order = ddata["leaves"]
    y_positions = [5 + i * 10 for i in range(len(leaf_order))]
    label_to_plin = dict(zip(labels, plin_codes))
    label_to_strain = dict(zip(labels, strain_clusters))

    for i, leaf_idx in enumerate(leaf_order):
        lbl = labels[leaf_idx]
        color = strain_cmap[label_to_strain[lbl]]
        ax_plin.text(0.05, y_positions[i], label_to_plin[lbl], fontsize=9,
                     fontfamily="monospace", fontweight="bold", color=color, va="center")

    legend_elements = [mpatches.Patch(color=strain_cmap[sc],
                       label=f"Strain {sc} (n={sum(1 for s in strain_clusters if s == sc)})")
                       for sc in unique_strains]
    ax_dendro.legend(handles=legend_elements, title="pLIN Strain (F)",
                     loc="upper right", fontsize=7, title_fontsize=8, framealpha=0.9)

    fig.suptitle(f"pLIN Cladogram (n={len(labels)})", fontsize=13, fontweight="bold", y=0.98)
    return fig


def plot_circular_cladogram(Z, labels, plin_codes, strain_clusters):
    strain_cmap = get_strain_cmap(strain_clusters)
    unique_strains = sorted(set(strain_clusters))
    n = len(labels)

    fig, ax = plt.subplots(figsize=(12, 12), subplot_kw={"polar": True})
    ddata = dendrogram(Z, labels=labels, no_plot=True)
    leaf_order = ddata["leaves"]

    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    angle_map = {leaf_order[i]: angles[i] for i in range(n)}

    max_dist = max(Z[:, 2]) * 1.15
    rs = 0.85

    node_pos = {}
    for i in range(n):
        node_pos[i] = (angle_map[i], 0.0)

    for idx, (c1, c2, dist, _) in enumerate(Z):
        c1, c2 = int(c1), int(c2)
        a1, d1 = node_pos[c1]
        a2, d2 = node_pos[c2]
        nid = n + idx
        avg_a = np.arctan2((np.sin(a1) + np.sin(a2)) / 2, (np.cos(a1) + np.cos(a2)) / 2)
        if avg_a < 0: avg_a += 2 * np.pi
        node_pos[nid] = (avg_a, dist)

        for a, d in [(a1, d1), (a2, d2)]:
            ax.plot([a, a], [(1 - d / max_dist) * rs, (1 - dist / max_dist) * rs],
                    color="#555", linewidth=1.5)
        a_min, a_max = min(a1, a2), max(a1, a2)
        arc = np.linspace(a_max, a_min + 2 * np.pi, 50) if (a_max - a_min) > np.pi else np.linspace(a_min, a_max, 50)
        ax.plot(arc, np.full_like(arc, (1 - dist / max_dist) * rs), color="#555", linewidth=1.5)

    for i in range(n):
        a = angle_map[i]
        sc = strain_clusters[i]
        color = strain_cmap[sc]
        ax.scatter(a, rs + 0.02, s=80, c=color, zorder=5, edgecolors="white", linewidth=0.5)
        rot = np.degrees(a) - 90
        ha = "left"
        if 90 < np.degrees(a) < 270:
            rot += 180
            ha = "right"
        ax.text(a, rs + 0.07, f"{labels[i]}  [{plin_codes[i]}]", fontsize=7,
                fontfamily="monospace", fontweight="bold", color=color,
                ha=ha, va="center", rotation=rot, rotation_mode="anchor")

    ax.set_ylim(0, rs + 0.22)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines["polar"].set_visible(False)
    ax.grid(False)
    ax.set_title(f"pLIN Circular Cladogram (n={n})", fontsize=13, fontweight="bold", pad=30)

    legend_elements = [mpatches.Patch(color=strain_cmap[sc], label=f"Strain {sc}") for sc in unique_strains]
    ax.legend(handles=legend_elements, title="pLIN Strain (F)", loc="lower left",
              bbox_to_anchor=(-0.05, -0.05), fontsize=7, title_fontsize=8, framealpha=0.9)
    return fig


def plot_cladogram_heatmap(Z, labels, plin_codes, cluster_assignments, records):
    fig = plt.figure(figsize=(16, max(8, len(labels) * 0.4)))
    ax_d = fig.add_axes([0.02, 0.08, 0.25, 0.82])
    ax_h = fig.add_axes([0.30, 0.08, 0.28, 0.82])
    ax_m = fig.add_axes([0.62, 0.08, 0.35, 0.82])

    ddata = dendrogram(Z, labels=labels, orientation="right", leaf_font_size=8,
                       ax=ax_d, color_threshold=0, above_threshold_color="#555")
    max_merge = Z[:, 2].max()
    ax_d.set_xlim(0, max_merge * 1.15)
    for bname, thresh in PLIN_THRESHOLDS.items():
        if thresh <= max_merge * 1.15:
            ax_d.axvline(x=thresh, color=THRESHOLD_COLORS[bname], linestyle="--", linewidth=1, alpha=0.6)
    ax_d.set_xlabel("Cosine Distance", fontsize=9)
    ax_d.set_title("Dendrogram", fontsize=10, fontweight="bold")

    leaf_order = ddata["leaves"]
    bin_labels = list(PLIN_THRESHOLDS.keys())
    heat = np.zeros((len(labels), len(bin_labels)))
    for j, b in enumerate(bin_labels):
        for i, li in enumerate(leaf_order):
            heat[i, j] = cluster_assignments[b][li]

    ax_h.imshow(heat, aspect="auto", cmap="tab20", interpolation="nearest")
    ax_h.set_xticks(range(len(bin_labels)))
    ax_h.set_xticklabels([f"{b}\n({PLIN_LEVEL_NAMES[b]})" for b in bin_labels], fontsize=8)
    ax_h.set_yticks(range(len(labels)))
    ax_h.set_yticklabels([labels[i] for i in leaf_order], fontsize=7)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            ax_h.text(j, i, str(int(heat[i, j])), ha="center", va="center", fontsize=7,
                      fontweight="bold", color="white",
                      path_effects=[pe.withStroke(linewidth=2, foreground="black")])
    ax_h.set_title("Cluster Assignments", fontsize=10, fontweight="bold")

    ax_m.axis("off")
    for j, h in enumerate(["Plasmid", "Length (bp)", "pLIN Code"]):
        ax_m.text(j * 0.35, len(labels) + 0.3, h, fontsize=8, fontweight="bold", ha="left")
    ax_m.axhline(y=len(labels), color="black", linewidth=0.8, xmin=0, xmax=0.95)
    for i, li in enumerate(leaf_order):
        y = len(labels) - 1 - i
        ax_m.text(0.0, y, records[li]["source_file"].replace(".fasta", ""), fontsize=7,
                  fontfamily="monospace", ha="left", va="center")
        ax_m.text(0.35, y, f"{records[li]['length']:,}", fontsize=7, ha="left", va="center")
        ax_m.text(0.70, y, plin_codes[li], fontsize=7, fontfamily="monospace",
                  fontweight="bold", color="#1565C0", ha="left", va="center")
    ax_m.set_xlim(-0.05, 1.1)
    ax_m.set_ylim(-1, len(labels) + 1)
    ax_m.set_title("Metadata", fontsize=10, fontweight="bold")

    fig.suptitle(f"pLIN Cladogram + Cluster Assignments (n={len(labels)})",
                 fontsize=13, fontweight="bold", y=0.97)
    return fig


def plot_cladogram_amr(Z, labels, plin_codes, strain_clusters, amr_df, records):
    """Cladogram with AMR gene presence/absence heatmap."""
    strain_cmap = get_strain_cmap(strain_clusters)
    unique_strains = sorted(set(strain_clusters))

    if len(amr_df) == 0:
        return plot_rectangular_cladogram(Z, labels, plin_codes, strain_clusters)

    amr_genes = sorted(amr_df[amr_df["Type"] == "AMR"]["Element symbol"].unique()) if "Type" in amr_df.columns else []
    stress_genes = sorted(amr_df[amr_df["Type"] == "STRESS"]["Element symbol"].unique()) if "Type" in amr_df.columns else []
    gene_list = amr_genes + stress_genes
    gene_types = ["AMR"] * len(amr_genes) + ["STRESS"] * len(stress_genes)

    if not gene_list:
        return plot_rectangular_cladogram(Z, labels, plin_codes, strain_clusters)

    n_genes = len(gene_list)
    n_plasmids = len(labels)
    label_to_source = {labels[i]: records[i]["source_file"].replace(".fasta", "").replace(".fa", "").replace(".fna", "")
                       for i in range(n_plasmids)}

    fig = plt.figure(figsize=(max(18, 12 + n_genes * 0.4), max(10, n_plasmids * 0.45)))
    ax_d = fig.add_axes([0.01, 0.10, 0.16, 0.78])
    ax_p = fig.add_axes([0.18, 0.10, 0.09, 0.78])
    ax_h = fig.add_axes([0.29, 0.10, min(0.50, n_genes * 0.025 + 0.1), 0.78])
    ax_l = fig.add_axes([0.82, 0.10, 0.16, 0.78])

    with plt.rc_context({"lines.linewidth": 2.0}):
        ddata = dendrogram(Z, labels=labels, orientation="right", leaf_font_size=8,
                           ax=ax_d, color_threshold=0, above_threshold_color="#444")
    max_merge = Z[:, 2].max()
    ax_d.set_xlim(0, max_merge * 1.15)
    for bname, thresh in PLIN_THRESHOLDS.items():
        if thresh <= max_merge * 1.15:
            ax_d.axvline(x=thresh, color=THRESHOLD_COLORS[bname], linestyle="--", linewidth=1, alpha=0.6)
    ax_d.set_xlabel("Cosine Distance", fontsize=8)
    ax_d.set_title("Dendrogram", fontsize=9, fontweight="bold")

    leaf_order = ddata["leaves"]
    leaf_colors = {labels[i]: strain_cmap[sc] for i, sc in enumerate(strain_clusters)}
    for lbl in ax_d.get_yticklabels():
        if lbl.get_text() in leaf_colors:
            lbl.set_color(leaf_colors[lbl.get_text()])
            lbl.set_fontweight("bold")

    # pLIN panel
    ax_p.set_ylim(ax_d.get_ylim())
    ax_p.set_xlim(0, 1)
    ax_p.axis("off")
    ax_p.set_title("pLIN", fontsize=9, fontweight="bold")
    y_pos = [5 + i * 10 for i in range(len(leaf_order))]
    label_to_plin = dict(zip(labels, plin_codes))
    label_to_strain = dict(zip(labels, strain_clusters))
    for i, li in enumerate(leaf_order):
        lbl = labels[li]
        ax_p.text(0.05, y_pos[i], label_to_plin[lbl], fontsize=7, fontfamily="monospace",
                  fontweight="bold", color=strain_cmap[label_to_strain[lbl]], va="center")

    # Heatmap
    heat = np.zeros((n_plasmids, n_genes))
    for i, li in enumerate(leaf_order):
        src = label_to_source[labels[li]]
        hits = amr_df[amr_df["source_file"] == src]
        for j, gene in enumerate(gene_list):
            if gene in hits["Element symbol"].values:
                heat[i, j] = 2 if gene_types[j] == "AMR" else 1

    cmap = ListedColormap(["#F5F5F5", "#FFE0B2", "#EF5350"])
    ax_h.imshow(heat, aspect="auto", cmap=cmap, interpolation="nearest", vmin=0, vmax=2)
    ax_h.set_xticks(range(n_genes))
    ax_h.set_xticklabels(gene_list, fontsize=6, rotation=65, ha="left", rotation_mode="anchor")
    for j, tick in enumerate(ax_h.get_xticklabels()):
        tick.set_color(TYPE_COLORS.get(gene_types[j], "#333"))
        tick.set_fontweight("bold")
    ax_h.set_yticks(range(n_plasmids))
    ax_h.set_yticklabels([labels[i] for i in leaf_order], fontsize=7)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            if heat[i, j] > 0:
                ax_h.text(j, i, "\u2713", ha="center", va="center", fontsize=6, fontweight="bold",
                          color="white" if heat[i, j] == 2 else "#E65100")
    ax_h.set_xticks(np.arange(-0.5, n_genes), minor=True)
    ax_h.set_yticks(np.arange(-0.5, n_plasmids), minor=True)
    ax_h.grid(which="minor", color="#E0E0E0", linewidth=0.5)
    ax_h.tick_params(which="minor", size=0)
    if amr_genes and stress_genes:
        ax_h.axvline(x=len(amr_genes) - 0.5, color="black", linewidth=2)
    ax_h.set_title("Gene Presence / Absence", fontsize=9, fontweight="bold")

    # Legend
    ax_l.axis("off")
    y = 0.95
    ax_l.text(0.05, y, "Strain Clusters", fontsize=8, fontweight="bold", transform=ax_l.transAxes)
    y -= 0.04
    for sc in unique_strains:
        n = sum(1 for s in strain_clusters if s == sc)
        ax_l.add_patch(mpatches.FancyBboxPatch((0.05, y - 0.01), 0.08, 0.025, transform=ax_l.transAxes,
                       boxstyle="round,pad=0.003", facecolor=strain_cmap[sc], edgecolor="none"))
        ax_l.text(0.16, y + 0.003, f"Strain {sc} (n={n})", fontsize=7,
                  transform=ax_l.transAxes, va="center")
        y -= 0.035

    fig.suptitle(f"pLIN Cladogram with AMR Profile (n={n_plasmids})",
                 fontsize=13, fontweight="bold", y=0.98)
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  SESSION STATE INITIALIZATION
# ══════════════════════════════════════════════════════════════════════════════

for key in ["records", "plin_df", "amr_df", "integrated_df", "Z", "labels",
            "plin_codes", "strain_clusters", "cluster_assignments", "analysis_done",
            "mobility_results", "outbreak_clusters", "active_thresholds",
            "linkage_method_used", "nt_results", "prodigal_genes_df", "prodigal_summary_df"]:
    if key not in st.session_state:
        st.session_state[key] = None
if "analysis_done" not in st.session_state:
    st.session_state.analysis_done = False

# Bacterial Buddy chat state
if "buddy_messages" not in st.session_state:
    st.session_state.buddy_messages = []
if "buddy_model" not in st.session_state:
    st.session_state.buddy_model = None


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN HEADER & UPLOAD
# ══════════════════════════════════════════════════════════════════════════════

st.title("🧬 pLIN Classifier")
st.caption("Plasmid Life Identification Number System — Upload FASTA files to begin")

# AMRFinderPlus detection (used in both upload and post-analysis views)
amr_binary, amr_db = detect_amrfinder()

# Prodigal detection
prodigal_binary = detect_prodigal()

# Ollama detection for Bacterial Buddy
ollama_available, ollama_models = detect_ollama()

if not st.session_state.analysis_done:
    st.divider()
    upload_col1, upload_col2 = st.columns([2, 1])
    with upload_col1:
        uploaded_files = st.file_uploader(
            "Upload plasmid FASTA files",
            type=["fasta", "fa", "fna"],
            accept_multiple_files=True,
            help="Upload one or more plasmid FASTA files (.fasta, .fa, .fna)",
            key="main_uploader",
        )
    with upload_col2:
        inc_type = st.selectbox(
            "Incompatibility Group",
            INC_GROUPS, index=0,
            help="'Auto-detect' uses a KNN classifier (96% accuracy) trained on 6,346 plasmids to identify Inc group per sequence",
        )
        linkage_method = st.selectbox(
            "Linkage Method",
            LINKAGE_METHODS, index=0,
            help="Single: traditional chaining (default). Complete: max distance, tighter clusters. Average: balanced. Weighted: WPGMA.",
        )
        use_adaptive = st.checkbox(
            "Adaptive thresholds",
            value=False,
            help="Calibrate pLIN thresholds per Inc group from training data distance distributions instead of fixed thresholds",
        )
        if amr_binary:
            run_amr = st.checkbox("Run AMRFinderPlus", value=True,
                                  help="Detect AMR, stress, and virulence genes using NCBI's AMRFinderPlus")
        else:
            st.warning("AMRFinderPlus not found", icon="⚠️")
            run_amr = False

        # Prodigal gene annotation (optional)
        if prodigal_binary:
            run_prodigal = st.checkbox(
                "Run Prodigal annotation",
                value=False,
                help="Annotate all genes using Prodigal (metagenomic mode). Shows gene count, coding density, and full gene table.",
            )
        else:
            run_prodigal = False
            st.info("Prodigal not found. Install: `conda install -c bioconda prodigal`", icon="🧬")

        # Nucleotide Transformer (optional LLM)
        if NT_AVAILABLE:
            use_nt = st.checkbox(
                "Use Nucleotide Transformer (LLM)",
                value=False,
                help="Use a genomic language model for enhanced Inc group and AMR class prediction. Requires trained probes (run train_nt_classifier.py first).",
            )
            if use_nt:
                nt_model_choice = st.selectbox(
                    "NT Model",
                    list(NT_MODELS.keys()), index=0,
                    help="Smaller models are faster; larger models may be more accurate.",
                )
                nt_device = detect_nt_device()
                st.caption(f"Device: {nt_device.upper()}")
            else:
                nt_model_choice = list(NT_MODELS.keys())[0]
                nt_device = "cpu"
        else:
            use_nt = False
            nt_model_choice = None
            nt_device = "cpu"
            st.info("Nucleotide Transformer not available. Install: `pip install transformers torch`", icon="🤖")

    if uploaded_files:
        st.info(f"📂 {len(uploaded_files)} file(s) uploaded — click **Run Analysis** below.")
        bcol1, bcol2, _ = st.columns([1, 1, 3])
        with bcol1:
            run_btn = st.button("▶ Run Analysis", type="primary", use_container_width=True)
        with bcol2:
            if st.button("🔄 Clear & Reset", use_container_width=True):
                for key in list(st.session_state.keys()):
                    del st.session_state[key]
                st.cache_data.clear()
                st.rerun()
    else:
        run_btn = False
        st.markdown("👆 **Upload FASTA files above** to get started.")
    st.divider()
else:
    # After analysis — show compact controls in main area
    uploaded_files = st.session_state.get("_uploaded_files", None)
    inc_type = st.session_state.get("_inc_type", INC_GROUPS[0])
    linkage_method = st.session_state.get("_linkage_method", "single")
    use_adaptive = st.session_state.get("_use_adaptive", False)
    use_nt = st.session_state.get("_use_nt", False)
    nt_model_choice = st.session_state.get("_nt_model_choice", None)
    nt_device = st.session_state.get("_nt_device", "cpu")
    run_amr = amr_binary is not None
    run_prodigal = st.session_state.get("_run_prodigal", False)
    run_btn = False


# ── Sidebar ───────────────────────────────────────────────────────────────────

with st.sidebar:
    st.title("🧬 pLIN")
    st.caption("Plasmid Life Identification Number")
    st.divider()

    if st.session_state.analysis_done:
        st.success("Analysis complete", icon="✅")
        df = st.session_state.plin_df
        if df is not None:
            st.metric("Plasmids", len(df))
            st.metric("Unique pLIN Codes", df["pLIN"].nunique())
            if "inc_confidence" in df.columns:
                inc_counts = df["inc_type"].value_counts()
                st.markdown("**Detected Inc Groups:**")
                for inc, cnt in inc_counts.items():
                    if inc == "Unknown/Novel":
                        st.markdown(f"- ⚠️ {inc}: {cnt}")
                    else:
                        st.markdown(f"- {inc}: {cnt}")
                # Show low-confidence warning in sidebar
                if "inc_is_low_confidence" in df.columns:
                    low_conf_count = df["inc_is_low_confidence"].sum()
                    if low_conf_count > 0:
                        st.warning(f"{low_conf_count} low-confidence")
            mob_df = st.session_state.get("mobility_results")
            if mob_df is not None and len(mob_df) > 0:
                mob_counts = mob_df["mobility"].value_counts().to_dict()
                conj = mob_counts.get("Conjugative", 0)
                if conj > 0:
                    st.markdown(f"**Mobility:** {conj} conjugative")
            outbreak = st.session_state.get("outbreak_clusters", [])
            if outbreak:
                st.markdown(f"**Outbreak clusters:** {len(outbreak)}")
            nt_res = st.session_state.get("nt_results")
            if nt_res and nt_res.get("inc_preds") is not None:
                st.markdown("**NT LLM:** enabled")
            prodigal_sum = st.session_state.get("prodigal_summary_df")
            if prodigal_sum is not None and len(prodigal_sum) > 0:
                total_genes = prodigal_sum["total_genes"].sum()
                st.markdown(f"**Prodigal:** {total_genes:,} genes")
        st.divider()

    if st.button("🔄 Clear & Reset", use_container_width=True, key="sidebar_reset"):
        for key in list(st.session_state.keys()):
            del st.session_state[key]
        st.cache_data.clear()
        st.rerun()

    if st.session_state.analysis_done:
        if st.button("📂 New Analysis", use_container_width=True, type="primary"):
            for key in list(st.session_state.keys()):
                del st.session_state[key]
            st.cache_data.clear()
            st.rerun()


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN ANALYSIS PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

if run_btn and uploaded_files:
    st.session_state._uploaded_files = uploaded_files
    st.session_state._inc_type = inc_type
    st.session_state._linkage_method = linkage_method
    st.session_state._use_adaptive = use_adaptive
    st.session_state._use_nt = use_nt
    st.session_state._nt_model_choice = nt_model_choice
    st.session_state._nt_device = nt_device
    st.session_state._run_prodigal = run_prodigal
    progress = st.progress(0, text="Starting analysis...")

    # Step 1: Parse sequences (+ Inc group auto-detection)
    progress.progress(5, text="Parsing FASTA files & detecting Inc groups..." if inc_type == "Auto-detect"
                      else "Parsing FASTA files...")
    try:
        records = parse_uploaded_fastas(uploaded_files, inc_type)
        if len(records) < 2:
            st.error("Need at least 2 sequences for clustering.")
            st.stop()
        st.session_state.records = records
    except Exception as e:
        st.error(f"Failed to parse FASTA files: {e}")
        st.stop()

    # Step 2: K-mer vectors
    progress.progress(15, text=f"Computing 4-mer vectors for {len(records)} plasmids...")
    sequences = [r["sequence"] for r in records]
    vectors = compute_kmer_vectors(tuple(sequences), k=4)

    # Step 2b: Adaptive threshold calibration (if enabled)
    active_thresholds = PLIN_THRESHOLDS
    if use_adaptive:
        progress.progress(25, text="Calibrating adaptive thresholds from training data...")
        calibrated = calibrate_inc_thresholds()
        if calibrated:
            # Use thresholds for the dominant Inc group in the upload
            inc_counts = Counter(r["inc_type"] for r in records)
            dominant_inc = inc_counts.most_common(1)[0][0]
            if dominant_inc in calibrated:
                active_thresholds = calibrated[dominant_inc]
                st.session_state.active_thresholds = active_thresholds
                st.session_state._calibrated_all = calibrated
                st.session_state._dominant_inc = dominant_inc

    # Step 3: Clustering & pLIN assignment
    progress.progress(40, text=f"Clustering ({linkage_method} linkage) & assigning pLIN codes...")
    plin_codes, cluster_assignments, Z, dist_condensed = assign_plin_codes(
        vectors, linkage_method=linkage_method, thresholds=active_thresholds
    )
    strain_clusters = list(cluster_assignments["F"])
    st.session_state.linkage_method_used = linkage_method

    # Step 4: Build results
    progress.progress(55, text="Building results table...")
    plin_df = build_results_df(records, plin_codes, cluster_assignments)
    labels = [r["source_file"].replace(".fasta", "").replace(".fa", "").replace(".fna", "")
              for r in records]

    st.session_state.plin_df = plin_df
    st.session_state.Z = Z
    st.session_state.labels = labels
    st.session_state.plin_codes = plin_codes
    st.session_state.strain_clusters = strain_clusters
    st.session_state.cluster_assignments = cluster_assignments

    # Step 4b: Nucleotide Transformer predictions (optional)
    if use_nt and NT_AVAILABLE:
        progress.progress(56, text="Loading Nucleotide Transformer model...")
        try:
            nt_model_id = NT_MODELS.get(nt_model_choice, list(NT_MODELS.values())[0])
            tokenizer, nt_model, nt_dev = load_nt_model(nt_model_id, nt_device)

            has_inc_probe = os.path.exists(NT_INC_PROBE_PATH)
            has_amr_probe = os.path.exists(NT_AMR_PROBE_PATH)

            if not has_inc_probe and not has_amr_probe:
                st.warning(
                    "NT probes not found. Run `python train_nt_classifier.py` first to train "
                    "the classifier probes from your training data.",
                    icon="🤖",
                )
                st.session_state.nt_results = None
            else:
                def nt_cb(pct):
                    progress.progress(int(56 + pct * 4),
                                      text=f"NT embeddings: {int(pct * 100)}%")

                nt_results = run_nt_predictions(sequences, tokenizer, nt_model, nt_dev, nt_cb)
                st.session_state.nt_results = nt_results
        except Exception as e:
            st.warning(f"Nucleotide Transformer failed: {e}. Continuing with KNN only.", icon="⚠️")
            st.session_state.nt_results = None

    # Step 5: AMRFinderPlus (optional)
    amr_df = pd.DataFrame()
    if run_amr and amr_binary:
        progress.progress(60, text="Running AMRFinderPlus...")

        def amr_cb(pct):
            progress.progress(int(60 + pct * 30), text=f"AMRFinderPlus: {int(pct * 100)}%")

        amr_df = run_amrfinder_on_files(uploaded_files, amr_binary, amr_db, amr_cb)

    st.session_state.amr_df = amr_df

    # Step 5b: Prodigal gene annotation (optional)
    prodigal_genes_df = pd.DataFrame()
    prodigal_summary_df = pd.DataFrame()
    if run_prodigal and prodigal_binary:
        progress.progress(75, text="Running Prodigal gene annotation...")

        def prodigal_cb(pct):
            progress.progress(int(75 + pct * 10), text=f"Prodigal: {int(pct * 100)}%")

        prodigal_genes_df, prodigal_summary_df = run_prodigal_on_files(
            uploaded_files, prodigal_binary, prodigal_cb
        )

    st.session_state.prodigal_genes_df = prodigal_genes_df
    st.session_state.prodigal_summary_df = prodigal_summary_df

    # Step 6: Integration
    progress.progress(90, text="Integrating results...")
    integrated_df = integrate_plin_amr(plin_df, amr_df)
    st.session_state.integrated_df = integrated_df

    # Step 7: Mobility prediction
    progress.progress(93, text="Predicting plasmid mobility...")
    mobility_results = []
    for rec in records:
        mob_class, mob_genes, mob_detail = classify_mobility(amr_df, rec["source_file"])
        mobility_results.append({
            "plasmid_id": rec["plasmid_id"],
            "mobility": mob_class,
            "mobility_genes": "; ".join(mob_genes) if mob_genes else "",
            "mobility_detail": mob_detail,
        })
    st.session_state.mobility_results = pd.DataFrame(mobility_results)

    # Step 8: Outbreak detection
    progress.progress(96, text="Scanning for outbreak clusters...")
    outbreak_clusters = detect_outbreak_clusters(plin_df, integrated_df)
    st.session_state.outbreak_clusters = outbreak_clusters

    progress.progress(100, text="Analysis complete!")
    st.session_state.analysis_done = True
    st.rerun()


# ══════════════════════════════════════════════════════════════════════════════
#  TABS
# ══════════════════════════════════════════════════════════════════════════════

tab_overview, tab_results, tab_clado, tab_amr, tab_epi, tab_buddy, tab_export = st.tabs(
    ["📋 Overview", "📊 Results", "🌳 Cladogram", "💊 AMR Analysis",
     "🔬 Epidemiology", "🦠 Bacterial Buddy", "📥 Export"]
)

# ── TAB 1: Overview ──────────────────────────────────────────────────────────

with tab_overview:
    st.header("pLIN — Plasmid Life Identification Number")
    st.markdown("""
    **pLIN** assigns each plasmid a six-position hierarchical code (`A.B.C.D.E.F`)
    based on tetranucleotide (4-mer) composition distances and single-linkage clustering.
    """)

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("pLIN Hierarchy")
        thresh_data = [{"Position": b, "Level": PLIN_LEVEL_NAMES[b],
                        "Threshold": f"d \u2264 {t:.3f}", "ANI": ANI_EQUIV[b]}
                       for b, t in PLIN_THRESHOLDS.items()]
        st.dataframe(pd.DataFrame(thresh_data), use_container_width=True, hide_index=True)

    with col2:
        st.subheader("How it works")
        st.markdown("""
        1. **Upload** plasmid FASTA sequences
        2. **Auto-detect** Inc group (KNN classifier, 96% accuracy)
        3. **Compute** tetranucleotide frequency vectors (256 features)
        4. **Calculate** pairwise cosine distances
        5. **Cluster** using single-linkage hierarchical method
        6. **Cut** tree at 6 thresholds → pLIN codes
        7. **Screen** for AMR genes with AMRFinderPlus *(optional)*
        """)

    if st.session_state.analysis_done:
        st.success("Analysis complete!")
        df = st.session_state.plin_df
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Plasmids", len(df))
        c2.metric("Unique pLIN Codes", df["pLIN"].nunique())
        c3.metric("Strain Clusters (F)", df["bin_F"].nunique())
        amr_count = st.session_state.integrated_df["AMR_count"].sum() if st.session_state.integrated_df is not None else 0
        c4.metric("AMR Detections", int(amr_count))

        # Show analysis parameters
        params_col1, params_col2 = st.columns(2)
        with params_col1:
            lm = st.session_state.get("linkage_method_used", "single")
            st.markdown(f"**Linkage method:** `{lm}`")
        with params_col2:
            at = st.session_state.get("active_thresholds")
            if at and at != PLIN_THRESHOLDS:
                dom_inc = st.session_state.get("_dominant_inc", "?")
                st.markdown(f"**Adaptive thresholds:** calibrated for **{dom_inc}**")
            else:
                st.markdown("**Thresholds:** default (fixed)")

        # Show adaptive threshold comparison if used
        if st.session_state.get("active_thresholds") and st.session_state.active_thresholds != PLIN_THRESHOLDS:
            with st.expander("Adaptive vs Default Thresholds"):
                at = st.session_state.active_thresholds
                comp_rows = []
                for level in "ABCDEF":
                    comp_rows.append({
                        "Level": f"{level} ({PLIN_LEVEL_NAMES[level]})",
                        "Default": f"{PLIN_THRESHOLDS[level]:.4f}",
                        "Adaptive": f"{at[level]:.6f}",
                        "Change": f"{((at[level] - PLIN_THRESHOLDS[level]) / PLIN_THRESHOLDS[level] * 100):+.1f}%",
                    })
                st.dataframe(pd.DataFrame(comp_rows), use_container_width=True, hide_index=True)

        # Inc group auto-detection summary
        if "inc_confidence" in df.columns:
            st.subheader("Inc Group Auto-Detection")

            # Check for low-confidence predictions (Unknown/Novel)
            has_low_conf = "inc_is_low_confidence" in df.columns
            if has_low_conf:
                low_conf_df = df[df["inc_is_low_confidence"] == True]
                low_conf_count = len(low_conf_df)
                if low_conf_count > 0:
                    st.error(
                        f"**{low_conf_count} plasmid(s) classified as Unknown/Novel Inc type** — "
                        f"confidence below {INC_CONFIDENCE_THRESHOLD*100:.0f}% threshold. "
                        f"These plasmids may represent novel Inc types not in the training data."
                    )
                    # Show low-confidence plasmids with their top 5 candidates
                    with st.expander(f"View Unknown/Novel plasmids ({low_conf_count})"):
                        low_conf_cols = ["plasmid_id", "inc_best_match", "inc_confidence"]
                        if "inc_top5_candidates" in df.columns:
                            low_conf_cols.append("inc_top5_candidates")
                        low_conf_display = low_conf_df[low_conf_cols].copy()
                        low_conf_display = low_conf_display.rename(columns={
                            "inc_best_match": "Best Match (Low Conf.)",
                            "inc_confidence": "Confidence",
                            "inc_top5_candidates": "Top 5 Candidates",
                        })
                        low_conf_display = low_conf_display.sort_values("Confidence", ascending=True)
                        st.dataframe(low_conf_display, use_container_width=True, hide_index=True)
                        st.caption(
                            "**Tip:** These plasmids did not match any known Inc group with sufficient confidence. "
                            "Consider: (1) manually verifying Inc type via replicon typing tools like PlasmidFinder, "
                            "(2) the plasmid may be a novel/rare Inc type, or (3) it may be a mosaic/hybrid plasmid."
                        )

            inc_summary = df.groupby("inc_type").agg(
                count=("plasmid_id", "count"),
                avg_confidence=("inc_confidence", "mean"),
                min_confidence=("inc_confidence", "min"),
            ).reset_index().rename(columns={
                "inc_type": "Inc Group",
                "count": "Plasmids",
                "avg_confidence": "Avg Confidence",
                "min_confidence": "Min Confidence",
            })
            inc_summary["Avg Confidence"] = inc_summary["Avg Confidence"].round(4)
            inc_summary["Min Confidence"] = inc_summary["Min Confidence"].round(4)
            st.dataframe(inc_summary, use_container_width=True, hide_index=True)

            # Warn about mixed Inc groups
            n_inc_groups = df["inc_type"].nunique()
            # Don't count Unknown/Novel as a real Inc group for the warning
            real_inc_groups = [g for g in df["inc_type"].unique() if g != "Unknown/Novel"]
            if len(real_inc_groups) > 1:
                st.info(
                    f"**{len(real_inc_groups)} Inc groups detected** in your upload. "
                    "pLIN clustering is computed across all plasmids together. "
                    "For within-group comparisons, upload files from a single Inc group."
                )

        # NT prediction comparison (if available)
        nt_res = st.session_state.get("nt_results")
        if nt_res and nt_res.get("inc_preds") is not None:
            st.subheader("Nucleotide Transformer LLM Predictions")
            nt_inc_preds = nt_res["inc_preds"]
            knn_preds = df["inc_type"].tolist()

            # Agreement rate
            agree = sum(1 for k, n in zip(knn_preds, nt_inc_preds) if k == n)
            total = len(knn_preds)
            agreement_pct = agree / total * 100 if total > 0 else 0

            nc1, nc2, nc3 = st.columns(3)
            nc1.metric("KNN vs NT Agreement", f"{agreement_pct:.1f}%")
            if nt_res.get("inc_cv_accuracy"):
                nc2.metric("NT Probe CV Accuracy", f"{nt_res['inc_cv_accuracy']:.1%}")
            nc3.metric("NT Embedding Dim", nt_res["embeddings"].shape[1])

            # Show disagreements
            disagreements = []
            for i, (k, n) in enumerate(zip(knn_preds, nt_inc_preds)):
                if k != n:
                    pid = df["plasmid_id"].iloc[i] if i < len(df) else f"seq_{i}"
                    disagreements.append({"Plasmid": pid, "KNN": k, "NT": n})
            if disagreements:
                with st.expander(f"KNN vs NT Disagreements ({len(disagreements)})"):
                    st.dataframe(pd.DataFrame(disagreements),
                                 use_container_width=True, hide_index=True)
            else:
                st.success("KNN and Nucleotide Transformer predictions agree on all plasmids.")

            # NT AMR predictions
            if nt_res.get("amr_preds") is not None:
                amr_classes = nt_res["amr_classes"]
                amr_preds = nt_res["amr_preds"]
                pos_counts = amr_preds.sum(axis=0)
                detected = [(c, int(n)) for c, n in zip(amr_classes, pos_counts) if n > 0]
                if detected:
                    st.markdown("**NT-predicted AMR drug classes:**")
                    amr_pred_df = pd.DataFrame(detected, columns=["Drug Class", "Plasmids"])
                    amr_pred_df = amr_pred_df.sort_values("Plasmids", ascending=False)
                    st.dataframe(amr_pred_df, use_container_width=True, hide_index=True)

        # Prodigal gene annotation results (if available)
        prodigal_summary = st.session_state.get("prodigal_summary_df")
        prodigal_genes = st.session_state.get("prodigal_genes_df")
        if prodigal_summary is not None and len(prodigal_summary) > 0:
            st.subheader("Prodigal Gene Annotation")

            # Summary metrics
            total_genes = prodigal_summary["total_genes"].sum()
            avg_coding_density = prodigal_summary["coding_density_pct"].mean()
            avg_gene_length = prodigal_summary["avg_gene_length_aa"].mean()

            pc1, pc2, pc3, pc4 = st.columns(4)
            pc1.metric("Total Genes", f"{total_genes:,}")
            pc2.metric("Avg Coding Density", f"{avg_coding_density:.1f}%")
            pc3.metric("Avg Gene Length", f"{avg_gene_length:.0f} aa")
            pc4.metric("Plasmids Annotated", len(prodigal_summary))

            # Per-plasmid summary table
            st.markdown("**Per-Plasmid Gene Statistics**")
            display_summary = prodigal_summary.rename(columns={
                "source_file": "Plasmid",
                "sequence_length": "Length (bp)",
                "total_genes": "Total Genes",
                "complete_genes": "Complete",
                "partial_genes": "Partial",
                "coding_density_pct": "Coding %",
                "avg_gene_length_aa": "Avg Gene (aa)",
            })
            st.dataframe(display_summary, use_container_width=True, hide_index=True)

            # Expandable full gene table
            if prodigal_genes is not None and len(prodigal_genes) > 0:
                with st.expander(f"View All Predicted Genes ({len(prodigal_genes)})"):
                    gene_display = prodigal_genes[["source_file", "gene_id", "start", "end",
                                                   "strand", "length_aa", "is_complete"]].copy()
                    gene_display = gene_display.rename(columns={
                        "source_file": "Plasmid",
                        "gene_id": "Gene ID",
                        "start": "Start",
                        "end": "End",
                        "strand": "Strand",
                        "length_aa": "Length (aa)",
                        "is_complete": "Complete",
                    })
                    st.dataframe(gene_display, use_container_width=True, hide_index=True)


# ── TAB 2: Results ───────────────────────────────────────────────────────────

with tab_results:
    if not st.session_state.analysis_done:
        st.info("Run analysis first to see results.")
    else:
        st.header("pLIN Assignments")

        df = st.session_state.integrated_df
        has_inc_conf = "inc_confidence" in st.session_state.plin_df.columns

        # Merge mobility data if available
        mob_df = st.session_state.get("mobility_results")
        if mob_df is not None and len(mob_df) > 0:
            if "mobility" not in df.columns:
                df = df.merge(mob_df[["plasmid_id", "mobility", "mobility_genes"]], on="plasmid_id", how="left")

        # Merge NT predictions if available
        nt_res = st.session_state.get("nt_results")
        if nt_res and nt_res.get("inc_preds") is not None and "nt_inc_type" not in df.columns:
            nt_inc_preds = nt_res["inc_preds"]
            if len(nt_inc_preds) == len(df):
                df["nt_inc_type"] = nt_inc_preds

        if df is not None and "AMR_genes" in df.columns:
            display_cols = ["plasmid_id", "inc_type"]
            if "nt_inc_type" in df.columns:
                display_cols.append("nt_inc_type")
            if has_inc_conf:
                display_cols.append("inc_confidence")
            display_cols += ["source_file", "length_bp", "pLIN"]
            if "mobility" in df.columns:
                display_cols.append("mobility")
            display_cols += ["AMR_count", "STRESS_count", "AMR_genes", "AMR_classes"]
            if "mobility_genes" in df.columns:
                display_cols.append("mobility_genes")
            # merge inc_confidence from plin_df if not already present
            if has_inc_conf and "inc_confidence" not in df.columns:
                df = df.merge(
                    st.session_state.plin_df[["plasmid_id", "inc_confidence"]],
                    on="plasmid_id", how="left",
                )
            display_cols = [c for c in display_cols if c in df.columns]
        else:
            df = st.session_state.plin_df
            if mob_df is not None and len(mob_df) > 0 and "mobility" not in df.columns:
                df = df.merge(mob_df[["plasmid_id", "mobility", "mobility_genes"]], on="plasmid_id", how="left")
            display_cols = ["plasmid_id", "inc_type"]
            if "nt_inc_type" in df.columns:
                display_cols.append("nt_inc_type")
            if has_inc_conf:
                display_cols.append("inc_confidence")
            display_cols += ["source_file", "length_bp", "pLIN"]
            if "mobility" in df.columns:
                display_cols.append("mobility")
            display_cols += ["bin_A", "bin_B", "bin_C", "bin_D", "bin_E", "bin_F"]

        # Search filter
        search = st.text_input("🔍 Search plasmid ID, pLIN code, or Inc group", "")
        if search:
            mask = (
                df["plasmid_id"].str.contains(search, case=False, na=False) |
                df["pLIN"].str.contains(search, case=False, na=False) |
                df["inc_type"].str.contains(search, case=False, na=False)
            )
            df_display = df[mask][display_cols]
        else:
            df_display = df[display_cols]

        st.dataframe(df_display, use_container_width=True, hide_index=True)

        # Summary stats
        with st.expander("Summary Statistics"):
            c1, c2 = st.columns(2)
            with c1:
                st.write("**pLIN Code Distribution**")
                st.dataframe(df["pLIN"].value_counts().reset_index().rename(
                    columns={"index": "pLIN", "pLIN": "pLIN Code", "count": "Count"}),
                    use_container_width=True, hide_index=True)
            with c2:
                st.write("**Sequence Length Statistics**")
                st.write(df["length_bp"].describe().to_frame("Value"))

            if has_inc_conf:
                st.write("**Per-Plasmid Inc Group Classification**")

                # Check for low-confidence predictions
                plin_df = st.session_state.plin_df
                has_low_conf = "inc_is_low_confidence" in plin_df.columns
                has_top5 = "inc_top5_candidates" in plin_df.columns

                # Build display columns
                inc_cols = ["plasmid_id", "inc_type", "inc_confidence"]
                if has_low_conf:
                    inc_cols.append("inc_is_low_confidence")
                if "inc_best_match" in plin_df.columns:
                    inc_cols.append("inc_best_match")
                if has_top5:
                    inc_cols.append("inc_top5_candidates")

                inc_detail = plin_df[inc_cols].copy()
                inc_detail = inc_detail.sort_values("inc_confidence", ascending=True)

                # Count and warn about low-confidence predictions
                if has_low_conf:
                    low_conf_count = inc_detail["inc_is_low_confidence"].sum()
                    if low_conf_count > 0:
                        st.warning(
                            f"**{low_conf_count} plasmid(s) have low-confidence Inc type predictions** "
                            f"(below {INC_CONFIDENCE_THRESHOLD*100:.0f}% threshold). "
                            f"These are labeled as 'Unknown/Novel' and may represent novel Inc types "
                            f"not present in the training data. See 'Top 5 Candidates' for the nearest matches."
                        )

                # Rename columns for better display
                rename_map = {
                    "inc_is_low_confidence": "Low Confidence?",
                    "inc_best_match": "Best Match",
                    "inc_top5_candidates": "Top 5 Candidates",
                }
                inc_detail = inc_detail.rename(columns=rename_map)

                st.dataframe(inc_detail, use_container_width=True, hide_index=True)


# ── TAB 3: Cladogram ────────────────────────────────────────────────────────

with tab_clado:
    if not st.session_state.analysis_done:
        st.info("Run analysis first to see cladograms.")
    else:
        st.header("Cladogram Visualization")

        viz_type = st.radio("Visualization type", ["Rectangular", "Circular", "Heatmap", "AMR Annotated"],
                            horizontal=True)

        Z = st.session_state.Z
        labels = st.session_state.labels
        plin_codes = st.session_state.plin_codes
        strain_clusters = st.session_state.strain_clusters

        with st.spinner("Generating cladogram..."):
            if viz_type == "Rectangular":
                fig = plot_rectangular_cladogram(Z, labels, plin_codes, strain_clusters)
            elif viz_type == "Circular":
                fig = plot_circular_cladogram(Z, labels, plin_codes, strain_clusters)
            elif viz_type == "Heatmap":
                fig = plot_cladogram_heatmap(Z, labels, plin_codes,
                                            st.session_state.cluster_assignments,
                                            st.session_state.records)
            else:  # AMR Annotated
                amr_df = st.session_state.amr_df
                if amr_df is not None and len(amr_df) > 0:
                    fig = plot_cladogram_amr(Z, labels, plin_codes, strain_clusters,
                                            amr_df, st.session_state.records)
                else:
                    st.warning("No AMR data available. Run analysis with AMRFinderPlus enabled.")
                    fig = plot_rectangular_cladogram(Z, labels, plin_codes, strain_clusters)

            st.pyplot(fig, use_container_width=True)

        # Download buttons
        col1, col2 = st.columns(2)
        with col1:
            st.download_button("📥 Download PNG", fig_to_bytes(fig, "png"),
                               f"pLIN_cladogram_{viz_type.lower()}.png", "image/png")
        with col2:
            st.download_button("📥 Download PDF", fig_to_bytes(fig, "pdf"),
                               f"pLIN_cladogram_{viz_type.lower()}.pdf", "application/pdf")
        plt.close(fig)


# ── TAB 4: AMR Analysis ─────────────────────────────────────────────────────

with tab_amr:
    if not st.session_state.analysis_done:
        st.info("Run analysis first to see AMR results.")
    elif st.session_state.amr_df is None or len(st.session_state.amr_df) == 0:
        st.warning("No AMR data. Run a new analysis with AMRFinderPlus enabled.")
    else:
        st.header("AMR Gene Analysis")
        amr_df = st.session_state.amr_df
        integrated = st.session_state.integrated_df

        # Metrics
        c1, c2, c3, c4 = st.columns(4)
        amr_only = amr_df[amr_df["Type"] == "AMR"] if "Type" in amr_df.columns else pd.DataFrame()
        stress_only = amr_df[amr_df["Type"] == "STRESS"] if "Type" in amr_df.columns else pd.DataFrame()
        c1.metric("Total Detections", len(amr_df))
        c2.metric("AMR Genes", len(amr_only))
        c3.metric("Stress Genes", len(stress_only))
        c4.metric("Plasmids with AMR", int((integrated["AMR_count"] > 0).sum()))

        # Interactive gene prevalence chart
        if len(amr_only) > 0:
            col1, col2 = st.columns(2)
            with col1:
                gene_freq = amr_only.groupby("Element symbol")["source_file"].nunique().sort_values(ascending=False)
                fig_bar = px.bar(
                    x=gene_freq.values, y=gene_freq.index,
                    orientation="h", title="AMR Gene Prevalence (# Plasmids)",
                    labels={"x": "Number of Plasmids", "y": "Gene"},
                    color_discrete_sequence=["#E53935"],
                )
                fig_bar.update_layout(height=400, yaxis=dict(autorange="reversed"))
                st.plotly_chart(fig_bar, use_container_width=True)

            with col2:
                if "Class" in amr_only.columns:
                    class_counts = amr_only["Class"].value_counts()
                    fig_pie = px.pie(
                        values=class_counts.values, names=class_counts.index,
                        title="AMR Drug Classes",
                        color_discrete_sequence=px.colors.qualitative.Set2,
                    )
                    fig_pie.update_layout(height=400)
                    st.plotly_chart(fig_pie, use_container_width=True)

        # Critical gene alerts
        if len(amr_only) > 0:
            critical_genes = {
                "Carbapenemases": ["blaKPC", "blaNDM", "blaOXA-48", "blaVIM", "blaIMP"],
                "ESBLs": ["blaCTX-M", "blaSHV-12", "blaTEM"],
                "Colistin resistance": ["mcr"],
            }
            alerts = []
            for category, prefixes in critical_genes.items():
                for prefix in prefixes:
                    matches = amr_only[amr_only["Element symbol"].str.startswith(prefix)]
                    if len(matches) > 0:
                        genes = matches["Element symbol"].unique()
                        n_plasmids = matches["source_file"].nunique()
                        alerts.append(f"**{category}**: {', '.join(genes)} detected in {n_plasmids} plasmid(s)")
            if alerts:
                st.warning("⚠️ **Critical Resistance Genes Detected**\n\n" + "\n\n".join(alerts))

        # Full AMR table
        with st.expander("Full AMR Detection Table"):
            st.dataframe(amr_df, use_container_width=True, hide_index=True)


# ── TAB 5: Epidemiology ──────────────────────────────────────────────────────

with tab_epi:
    if not st.session_state.analysis_done:
        st.info("Run analysis first to see epidemiological insights.")
    else:
        st.header("Epidemiological Analysis")

        epi_col1, epi_col2 = st.columns(2)

        # ─── Mobility Prediction ───
        with epi_col1:
            st.subheader("Plasmid Mobility Prediction")
            mob_df = st.session_state.get("mobility_results")
            if mob_df is not None and len(mob_df) > 0:
                mob_counts = mob_df["mobility"].value_counts()
                # Color-coded metrics
                mob_colors = {"Conjugative": "#E53935", "Mobilizable": "#FB8C00",
                              "Non-mobilizable": "#43A047", "Unknown": "#9E9E9E"}
                mc1, mc2, mc3 = st.columns(3)
                mc1.metric("Conjugative", int(mob_counts.get("Conjugative", 0)))
                mc2.metric("Mobilizable", int(mob_counts.get("Mobilizable", 0)))
                mc3.metric("Non-mobilizable", int(mob_counts.get("Non-mobilizable", 0)))

                # Pie chart
                fig_mob = px.pie(
                    values=mob_counts.values, names=mob_counts.index,
                    title="Mobility Classification",
                    color=mob_counts.index,
                    color_discrete_map=mob_colors,
                )
                fig_mob.update_layout(height=350)
                st.plotly_chart(fig_mob, use_container_width=True)

                # Detail table
                with st.expander("Mobility Details"):
                    st.dataframe(mob_df, use_container_width=True, hide_index=True)

                if int(mob_counts.get("Conjugative", 0)) > 0:
                    conj_plasmids = mob_df[mob_df["mobility"] == "Conjugative"]["plasmid_id"].tolist()
                    st.warning(
                        f"**{len(conj_plasmids)} conjugative plasmid(s) detected** — "
                        "these can self-transfer to other bacteria via conjugation, "
                        "posing higher risk for AMR dissemination."
                    )
            else:
                st.info("No mobility data. Run analysis with AMRFinderPlus enabled.")

        # ─── Outbreak Detection ───
        with epi_col2:
            st.subheader("Outbreak / Clone Detection")
            outbreak_clusters = st.session_state.get("outbreak_clusters", [])

            if outbreak_clusters:
                st.metric("Suspected Outbreak Clusters", len(outbreak_clusters))

                for i, cluster in enumerate(outbreak_clusters):
                    risk_color = "🔴" if cluster["risk_level"] == "HIGH" else "🟡"
                    with st.expander(
                        f"{risk_color} Cluster {i+1}: pLIN {cluster['pLIN']} "
                        f"({cluster['n_plasmids']} plasmids, {cluster['n_amr_genes']} AMR genes)"
                    ):
                        st.markdown(f"**Risk level:** {cluster['risk_level']}")
                        st.markdown(f"**Strain cluster (F):** {cluster['strain_cluster']}")
                        st.markdown(f"**Shared AMR genes:** {', '.join(cluster['amr_genes'])}")
                        st.markdown(f"**Plasmids:**")
                        for p in cluster["plasmids"]:
                            st.markdown(f"  - `{p}`")

                if any(c["risk_level"] == "HIGH" for c in outbreak_clusters):
                    st.error(
                        "**HIGH-RISK outbreak cluster(s) detected.** "
                        "Plasmids sharing identical strain-level pLIN codes AND "
                        "the same AMR resistance profile may indicate clonal spread."
                    )
            else:
                st.info(
                    "No outbreak clusters detected. Outbreak detection flags groups of "
                    "plasmids sharing the same pLIN strain code (F-level) AND identical "
                    "AMR resistance profiles."
                )

        # ─── Risk Summary ───
        st.divider()
        st.subheader("Dissemination Risk Summary")

        integrated = st.session_state.integrated_df
        mob_df = st.session_state.get("mobility_results")

        if integrated is not None and mob_df is not None and len(mob_df) > 0:
            risk_df = integrated[["plasmid_id", "pLIN", "inc_type", "AMR_count"]].merge(
                mob_df[["plasmid_id", "mobility"]], on="plasmid_id", how="left"
            )

            # High risk: conjugative + AMR genes
            high_risk = risk_df[(risk_df["mobility"] == "Conjugative") & (risk_df["AMR_count"] > 0)]
            moderate_risk = risk_df[(risk_df["mobility"] == "Mobilizable") & (risk_df["AMR_count"] > 0)]

            rc1, rc2, rc3 = st.columns(3)
            rc1.metric("High Risk", len(high_risk),
                       help="Conjugative plasmids carrying AMR genes")
            rc2.metric("Moderate Risk", len(moderate_risk),
                       help="Mobilizable plasmids carrying AMR genes")
            rc3.metric("Lower Risk", len(risk_df) - len(high_risk) - len(moderate_risk),
                       help="Non-mobilizable or no AMR genes")

            if len(high_risk) > 0:
                st.warning(
                    f"**{len(high_risk)} high-risk plasmid(s):** conjugative AND carrying AMR genes. "
                    "These represent the highest priority for infection control surveillance."
                )
                st.dataframe(high_risk[["plasmid_id", "pLIN", "inc_type", "AMR_count", "mobility"]],
                             use_container_width=True, hide_index=True)


# ── TAB 6: Bacterial Buddy ───────────────────────────────────────────────────

with tab_buddy:
    st.header("🦠 Bacterial Buddy")
    st.caption("Your AI assistant for plasmid biology and AMR analysis")

    if not ollama_available:
        st.warning(
            "**Ollama not detected.** Bacterial Buddy requires Ollama running locally.\n\n"
            "**Setup instructions:**\n"
            "1. Install Ollama: https://ollama.ai\n"
            "2. Start Ollama: `ollama serve`\n"
            "3. Pull a model: `ollama pull llama3.2`\n"
            "4. Refresh this page",
            icon="🤖"
        )
        st.info(
            "**Why Ollama?** It runs LLMs locally on your machine — no API keys, no data leaves your computer, "
            "and it's completely free. Perfect for sensitive research data!",
            icon="💡"
        )
    else:
        # Model selection
        col1, col2 = st.columns([2, 1])
        with col1:
            # Filter to models that are actually installed
            available_models = [m for m in OLLAMA_MODELS if m in ollama_models] or ollama_models[:5]
            if available_models:
                selected_model = st.selectbox(
                    "Select AI Model",
                    available_models,
                    index=0,
                    help="Smaller models (llama3.2, phi3) are faster. Larger models (mixtral) may give better answers.",
                )
                st.session_state.buddy_model = selected_model
            else:
                st.warning("No models found. Run `ollama pull llama3.2` to download a model.")
                selected_model = None

        with col2:
            if st.button("🗑️ Clear Chat", use_container_width=True):
                st.session_state.buddy_messages = []
                st.rerun()

        # Show analysis status
        if st.session_state.analysis_done:
            st.success("Analysis data loaded. Bacterial Buddy can answer questions about your results!", icon="✅")
        else:
            st.info("Run an analysis first for context-aware answers, or ask general plasmid biology questions.", icon="💡")

        st.divider()

        # Chat container
        chat_container = st.container()

        # Display chat history
        with chat_container:
            for message in st.session_state.buddy_messages:
                if message["role"] == "user":
                    with st.chat_message("user", avatar="👤"):
                        st.markdown(message["content"])
                else:
                    with st.chat_message("assistant", avatar="🦠"):
                        st.markdown(message["content"])

        # Chat input
        if selected_model:
            if prompt := st.chat_input("Ask Bacterial Buddy about your plasmids..."):
                # Add user message
                st.session_state.buddy_messages.append({"role": "user", "content": prompt})

                # Display user message
                with chat_container:
                    with st.chat_message("user", avatar="👤"):
                        st.markdown(prompt)

                # Build context from analysis
                context = build_analysis_context(st.session_state)

                # Get and display assistant response
                with chat_container:
                    with st.chat_message("assistant", avatar="🦠"):
                        response_placeholder = st.empty()
                        full_response = ""

                        # Stream the response
                        for chunk in stream_chat_with_ollama(
                            st.session_state.buddy_messages,
                            selected_model,
                            context
                        ):
                            full_response += chunk
                            response_placeholder.markdown(full_response + "▌")

                        response_placeholder.markdown(full_response)

                # Save assistant message
                st.session_state.buddy_messages.append({"role": "assistant", "content": full_response})

        # Suggested questions
        st.divider()
        st.markdown("**💡 Try asking:**")
        suggestions = [
            "What Inc groups are in my data and what do they mean?",
            "Which plasmids are most likely to spread resistance?",
            "Explain the AMR genes found in my analysis",
            "What is the clinical significance of conjugative plasmids?",
            "How does pLIN classification work?",
        ]
        if st.session_state.analysis_done:
            suggestions = [
                "Summarize my analysis results",
                "Which plasmids should I be most concerned about?",
                "Are there any potential outbreak clusters?",
                "Explain the Inc groups detected in my samples",
                "What AMR genes are most prevalent and why does it matter?",
            ]

        cols = st.columns(2)
        for i, suggestion in enumerate(suggestions):
            with cols[i % 2]:
                if st.button(f"💬 {suggestion}", key=f"suggest_{i}", use_container_width=True):
                    # Trigger the question
                    st.session_state.buddy_messages.append({"role": "user", "content": suggestion})
                    st.rerun()


# ── TAB 7: Export ────────────────────────────────────────────────────────────

with tab_export:
    if not st.session_state.analysis_done:
        st.info("Run analysis first to export results.")
    else:
        st.header("Export Results")

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Tables")
            # pLIN results
            plin_csv = st.session_state.plin_df.to_csv(sep="\t", index=False).encode()
            st.download_button("📥 pLIN Assignments (TSV)", plin_csv,
                               "pLIN_assignments.tsv", "text/tab-separated-values")

            # Integrated results
            if st.session_state.integrated_df is not None:
                int_csv = st.session_state.integrated_df.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Integrated pLIN + AMR (TSV)", int_csv,
                                   "pLIN_AMR_integrated.tsv", "text/tab-separated-values")

            # Raw AMR
            if st.session_state.amr_df is not None and len(st.session_state.amr_df) > 0:
                amr_csv = st.session_state.amr_df.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 AMR Detections (TSV)", amr_csv,
                                   "amrfinder_results.tsv", "text/tab-separated-values")

            # Mobility results
            mob_df = st.session_state.get("mobility_results")
            if mob_df is not None and len(mob_df) > 0:
                mob_csv = mob_df.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Mobility Predictions (TSV)", mob_csv,
                                   "mobility_predictions.tsv", "text/tab-separated-values")

            # Prodigal results
            prodigal_summary = st.session_state.get("prodigal_summary_df")
            prodigal_genes = st.session_state.get("prodigal_genes_df")
            if prodigal_summary is not None and len(prodigal_summary) > 0:
                prodigal_sum_csv = prodigal_summary.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Prodigal Summary (TSV)", prodigal_sum_csv,
                                   "prodigal_summary.tsv", "text/tab-separated-values")
            if prodigal_genes is not None and len(prodigal_genes) > 0:
                prodigal_genes_csv = prodigal_genes.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Prodigal Genes (TSV)", prodigal_genes_csv,
                                   "prodigal_genes.tsv", "text/tab-separated-values")

        with col2:
            st.subheader("Figures")
            Z = st.session_state.Z
            labels = st.session_state.labels
            plin_codes = st.session_state.plin_codes
            strain_clusters = st.session_state.strain_clusters

            for name, func, args in [
                ("Rectangular Cladogram", plot_rectangular_cladogram,
                 (Z, labels, plin_codes, strain_clusters)),
                ("Circular Cladogram", plot_circular_cladogram,
                 (Z, labels, plin_codes, strain_clusters)),
            ]:
                fig = func(*args)
                st.download_button(f"📥 {name} (PNG)", fig_to_bytes(fig),
                                   f"{name.lower().replace(' ', '_')}.png", "image/png")
                plt.close(fig)

            # AMR cladogram if available
            if st.session_state.amr_df is not None and len(st.session_state.amr_df) > 0:
                fig = plot_cladogram_amr(Z, labels, plin_codes, strain_clusters,
                                        st.session_state.amr_df, st.session_state.records)
                st.download_button("📥 AMR Cladogram (PNG)", fig_to_bytes(fig),
                                   "cladogram_amr.png", "image/png")
                plt.close(fig)

        st.divider()

        # ZIP bundle
        st.subheader("Download All (ZIP)")
        if st.button("📦 Generate ZIP Bundle"):
            with st.spinner("Creating ZIP..."):
                zip_buf = io.BytesIO()
                with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
                    zf.writestr("pLIN_assignments.tsv",
                                st.session_state.plin_df.to_csv(sep="\t", index=False))
                    if st.session_state.integrated_df is not None:
                        zf.writestr("pLIN_AMR_integrated.tsv",
                                    st.session_state.integrated_df.to_csv(sep="\t", index=False))
                    if st.session_state.amr_df is not None and len(st.session_state.amr_df) > 0:
                        zf.writestr("amrfinder_results.tsv",
                                    st.session_state.amr_df.to_csv(sep="\t", index=False))
                    mob_df = st.session_state.get("mobility_results")
                    if mob_df is not None and len(mob_df) > 0:
                        zf.writestr("mobility_predictions.tsv",
                                    mob_df.to_csv(sep="\t", index=False))
                    outbreak = st.session_state.get("outbreak_clusters", [])
                    if outbreak:
                        import json
                        zf.writestr("outbreak_clusters.json",
                                    json.dumps(outbreak, indent=2))
                    # NT predictions
                    nt_res = st.session_state.get("nt_results")
                    if nt_res and nt_res.get("inc_preds") is not None:
                        nt_df = pd.DataFrame({
                            "plasmid_id": st.session_state.plin_df["plasmid_id"],
                            "knn_inc_type": st.session_state.plin_df["inc_type"],
                            "nt_inc_type": nt_res["inc_preds"],
                        })
                        zf.writestr("nt_predictions.tsv",
                                    nt_df.to_csv(sep="\t", index=False))
                    # Prodigal results
                    prodigal_summary = st.session_state.get("prodigal_summary_df")
                    prodigal_genes = st.session_state.get("prodigal_genes_df")
                    if prodigal_summary is not None and len(prodigal_summary) > 0:
                        zf.writestr("prodigal_summary.tsv",
                                    prodigal_summary.to_csv(sep="\t", index=False))
                    if prodigal_genes is not None and len(prodigal_genes) > 0:
                        zf.writestr("prodigal_genes.tsv",
                                    prodigal_genes.to_csv(sep="\t", index=False))
                    for name, func, args in [
                        ("cladogram_rectangular", plot_rectangular_cladogram,
                         (Z, labels, plin_codes, strain_clusters)),
                        ("cladogram_circular", plot_circular_cladogram,
                         (Z, labels, plin_codes, strain_clusters)),
                    ]:
                        fig = func(*args)
                        zf.writestr(f"{name}.png", fig_to_bytes(fig))
                        zf.writestr(f"{name}.pdf", fig_to_bytes(fig, "pdf"))
                        plt.close(fig)

                    # AMR cladogram in ZIP
                    if st.session_state.amr_df is not None and len(st.session_state.amr_df) > 0:
                        fig = plot_cladogram_amr(Z, labels, plin_codes, strain_clusters,
                                                st.session_state.amr_df, st.session_state.records)
                        zf.writestr("cladogram_amr.png", fig_to_bytes(fig))
                        zf.writestr("cladogram_amr.pdf", fig_to_bytes(fig, "pdf"))
                        plt.close(fig)

                zip_buf.seek(0)
                st.download_button("📥 Download ZIP", zip_buf.getvalue(),
                                   "pLIN_results.zip", "application/zip")
