#!/usr/bin/env python3
"""
pLIN Classifier — Streamlit GUI Application
Plasmid Lineage Identification Number system with AMRFinderPlus integration.
Run: streamlit run plin_app.py

Copyright (C) 2025 Basil Xavier Britto
Licensed under GPL-3.0 with mandatory citation clause.
See LICENSE and CITATION.cff for details.

CITATION REQUIRED: Any use of this software in publications or derivative
works must cite:
    Xavier, B. (2025). pLIN: A Plasmid Lineage Identification Number System
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
import shutil
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
from scipy.spatial.distance import pdist, squareform, cdist
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from Bio import SeqIO
from sklearn.neighbors import KNeighborsClassifier
from collections import Counter

# ── Page Configuration ────────────────────────────────────────────────────────

_LOGO_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")
_FAVICON = os.path.join(_LOGO_DIR, "pLIN_favicon.png")

st.set_page_config(
    page_title="Generating pLIN number — a digital ID for your plasmid",
    page_icon=_FAVICON if os.path.exists(_FAVICON) else "🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Constants ─────────────────────────────────────────────────────────────────

PLIN_THRESHOLDS = {
    "A": 0.150, "B": 0.100, "C": 0.050,
    "D": 0.020, "E": 0.010, "F": 0.001,
}

PLIN_LEVEL_NAMES = {
    "A": "L1", "B": "L2", "C": "L3",
    "D": "L4", "E": "L5", "F": "L6",
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

# Paths to precomputed Inc-group classifier data
_APP_DIR = os.path.dirname(os.path.abspath(__file__))
CLASSIFIER_PATH = os.path.join(_APP_DIR, "data", "inc_classifier.npz")
CENTROID_PATH = os.path.join(_APP_DIR, "data", "inc_centroids.npz")

def _get_inc_groups():
    """Load Inc group names from classifier data, with Auto-detect and Other."""
    groups = ["Auto-detect"]
    if os.path.exists(CLASSIFIER_PATH):
        data = np.load(CLASSIFIER_PATH, allow_pickle=True)
        groups.extend(sorted(str(g) for g in data["group_names"]))
    else:
        groups.extend(["IncFII", "IncN", "IncX1"])  # fallback
    groups.append("Other")
    return groups

INC_GROUPS = _get_inc_groups()

LINKAGE_METHODS = ["single", "complete", "average", "weighted"]

# Confidence threshold for Inc group classification
# Below this threshold, plasmids are flagged as "Unknown/Novel" Inc type
INC_CONFIDENCE_THRESHOLD = 0.40  # 40% confidence minimum

# Threshold for detecting multiple Inc types (multi-replicon plasmids)
# If 2+ Inc types have confidence >= this threshold, flag as "Multiple Inc"
# Lowered from 0.25 → 0.15: large multi-replicon plasmids (e.g. IncHI2/IncN co-carrying
# blaVIM-1) split KNN probability ~80/20 between the two replicon types; the 25% threshold
# missed the secondary replicon for 4 of 8 Swiss VIM-1 isolates.
MULTI_INC_THRESHOLD = 0.15  # 15% minimum for secondary Inc types

# PlasmidFinder replicon BLAST — used as a secondary signal for large plasmids
# where KNN gives ≥95% to a single group (all 5 neighbours same group → misses secondary)
PLASMIDFINDER_DB_DIR = os.path.expanduser("~/plasmidfinder_db")
REPLICON_BLAST_SIZE_THRESHOLD = 100_000   # only run for sequences > 100 kb
REPLICON_BLAST_IDENTITY = 80.0            # min % identity for replicon hit
REPLICON_BLAST_COVERAGE = 60.0           # min query coverage % for replicon hit

# Minimum sequence length for reliable 4-mer classification
# Plasmids shorter than this have high stochastic variance in k-mer profiles
SHORT_PLASMID_THRESHOLD = 5000  # 5 kb — warn users about unreliable pLIN codes
CHROMOSOMAL_THRESHOLD = 500000  # 500 kb — likely chromosomal, not plasmid

# Mobility/conjugation marker genes detectable from AMRFinderPlus output
MOBILITY_GENES = {
    "conjugative": {
        "tra": "Transfer (F-type conjugation)",
        "trb": "Transfer (IncP/Ti-type T4SS)",
        "trw": "Transfer (IncW-type T4SS)",
        "virB": "Type IV secretion system (T4SS)",
        "virD": "Relaxase/coupling protein (T4SS)",
        "pil": "Type IV pilus (conjugation)",
    },
    "mobilizable": {
        "mob": "Mobilization protein",
        "oriT": "Origin of transfer",
        "nic/nik": "Nickase (relaxase)",
        "MOBF/MOBH/MOBP/MOBQ/MOBC/MOBV": "Relaxase MOB families",
    },
}

# All mobility gene prefixes for scanning
MOBILITY_PREFIXES_CONJUGATIVE = [
    # tra genes (F-type conjugation system)
    "traA", "traB", "traC", "traD", "traE", "traF",
    "traG", "traH", "traI", "traJ", "traK", "traL",
    "traM", "traN", "traO", "traP", "traQ", "traR",
    "traS", "traT", "traU", "traV", "traW", "traX", "traY",
    # trb genes (IncP-type / Ti-type T4SS)
    "trbA", "trbB", "trbC", "trbD", "trbE", "trbF",
    "trbG", "trbH", "trbI", "trbJ", "trbK", "trbL", "trbM", "trbN",
    # trw genes (IncW-type T4SS)
    "trwA", "trwB", "trwC", "trwD", "trwE", "trwF", "trwG", "trwH",
    "trwI", "trwJ", "trwK", "trwL", "trwM", "trwN",
    # virB/virD genes (Agrobacterium-like T4SS, also in conjugative plasmids)
    "virB1", "virB2", "virB3", "virB4", "virB5", "virB6", "virB7",
    "virB8", "virB9", "virB10", "virB11", "virD2", "virD4",
    # Type IV pilus
    "pilX",
    # Transfer accessory
    "taxC",
]
MOBILITY_PREFIXES_MOBILIZABLE = [
    # mob genes (mobilization proteins)
    "mobA", "mobB", "mobC", "mobD", "mobE", "mobF",
    # nik genes (nickase/relaxase)
    "nikA", "nikB", "nikC", "nikD", "nikE",
    # Relaxase MOB families (as classified by MOBscan)
    "MOBF", "MOBH", "MOBP", "MOBQ", "MOBC", "MOBV",
    # oriT-associated
    "oriT",
    # Mobilization-associated primase/nickase
    "mps",
]

# Keyword patterns for scanning gene names/descriptions
MOBILITY_KEYWORDS_CONJUGATIVE = [
    "conjugal transfer", "conjugative transfer", "conjugation",
    "type iv secretion", "type 4 secretion", "t4ss",
    "mating pair formation",
    "coupling protein", "t4cp",
    "dna transfer", "sex pilus",
]
MOBILITY_KEYWORDS_MOBILIZABLE = [
    "mobilization", "mobilisation",
    "relaxase", "relaxosome",
    "nickase", "origin of transfer", "orit",
    "mob family",
]

# ── CRISPR Host Inference Constants ──────────────────────────────────────────

CRISPR_BLAST_IDENTITY_MIN = 95.0      # Minimum % identity for spacer-plasmid hit
CRISPR_BLAST_ALIGNMENT_MIN = 25       # Minimum alignment length (bp)
CRISPR_BLAST_MISMATCH_MAX = 1         # Maximum allowed mismatches
CRISPR_BLAST_GAPS_MAX = 0             # No gaps allowed
CRISPR_BLAST_EVALUE = 1e-5            # E-value cutoff for blastn-short
CRISPR_SOFTMAX_TEMPERATURE = 1.0      # Softmax temperature (higher = flatter distribution)

REFERENCE_FASTA_PATH = os.path.join(_APP_DIR, "sequences.fasta")
REFERENCE_BLASTDB_PATH = os.path.join(_APP_DIR, "reference", "plasmid_blastdb")

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
    """Load precomputed Inc-group KNN classifier (91.1% CV accuracy, 28 groups)."""
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


@st.cache_resource(show_spinner=False)
def load_cv_metrics():
    """Load per-Inc-group CV metrics from classifier data."""
    import json as _json
    if os.path.exists(CLASSIFIER_PATH):
        data = np.load(CLASSIFIER_PATH, allow_pickle=True)
        cv_metrics = {}
        if "cv_metrics" in data:
            raw = str(data["cv_metrics"][0])
            cv_metrics = _json.loads(raw)
        cv_accuracy = float(data["cv_accuracy"][0]) if "cv_accuracy" in data else 0.0
        confusion_pairs = []
        if "confusion_pairs" in data:
            for cp in data["confusion_pairs"]:
                s = str(cp)
                if "|" in s:
                    parts = s.split("|")
                    confusion_pairs.append({
                        "inc_a": parts[0], "inc_b": parts[1],
                        "errors": int(parts[2]),
                        "pct_a": float(parts[3]), "pct_b": float(parts[4]),
                    })
        return {"per_class": cv_metrics, "accuracy": cv_accuracy,
                "confusion_pairs": confusion_pairs}
    return {"per_class": {}, "accuracy": 0.0, "confusion_pairs": []}


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


# ── Input Quality Validation ─────────────────────────────────────────────────

def validate_sequence_quality(sequence, plasmid_id):
    """Check input sequence quality before classification.

    Returns list of dicts with keys: level, check, message, value.
    Levels: 'error' (unreliable), 'warning' (interpret with caution), 'info'.
    """
    warnings = []
    seq_upper = sequence.upper()
    seq_len = len(seq_upper)
    if seq_len == 0:
        return [{"level": "error", "check": "empty",
                 "message": f"{plasmid_id}: empty sequence", "value": 0}]

    # 1. N-content check
    n_count = seq_upper.count("N")
    n_pct = 100.0 * n_count / seq_len
    if n_pct > 20:
        warnings.append({"level": "error", "check": "N-content",
                          "message": f"{plasmid_id}: {n_pct:.1f}% N bases — assembly too fragmented for reliable 4-mer profiling",
                          "value": round(n_pct, 1)})
    elif n_pct > 5:
        warnings.append({"level": "warning", "check": "N-content",
                          "message": f"{plasmid_id}: {n_pct:.1f}% N bases — may reduce 4-mer accuracy",
                          "value": round(n_pct, 1)})

    # 2. GC-content check (25–70% covers Enterobacterales + Gram-positive plasmids)
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    valid_bases = sum(seq_upper.count(b) for b in "ACGT")
    gc_pct = 100.0 * gc_count / valid_bases if valid_bases > 0 else 0
    if gc_pct < 25 or gc_pct > 70:
        warnings.append({"level": "warning", "check": "GC-content",
                          "message": f"{plasmid_id}: GC = {gc_pct:.1f}% (expected 25–70% for bacterial plasmids)",
                          "value": round(gc_pct, 1)})

    # 3. Non-ACGTN characters
    non_standard = sum(1 for c in seq_upper if c not in "ACGTN")
    if non_standard > 0:
        non_pct = 100.0 * non_standard / seq_len
        warnings.append({"level": "warning", "check": "non-ACGTN",
                          "message": f"{plasmid_id}: {non_standard} non-ACGTN characters ({non_pct:.2f}%)",
                          "value": non_standard})

    # 4. Low-complexity check (Shannon entropy of 4-mer distribution)
    if seq_len >= 100:
        vec = _kmer_vector_single(sequence)
        nonzero = vec[vec > 0]
        if len(nonzero) > 0:
            entropy = -np.sum(nonzero * np.log2(nonzero))
            # Max entropy for 256 bins = log2(256) = 8.0; typical plasmid > 6.0
            if entropy < 3.5:
                warnings.append({"level": "warning", "check": "low-complexity",
                                  "message": f"{plasmid_id}: 4-mer entropy = {entropy:.2f} bits (very low — repetitive sequence)",
                                  "value": round(entropy, 2)})
            elif entropy < 5.0:
                warnings.append({"level": "info", "check": "low-complexity",
                                  "message": f"{plasmid_id}: 4-mer entropy = {entropy:.2f} bits (below average)",
                                  "value": round(entropy, 2)})

    return warnings


# ── Plasmid vs Chromosome Contig Classification ─────────────────────────────

def classify_contigs_plasmid_vs_chromosome(records, vectors=None):
    """Classify each contig as 'plasmid' or 'chromosome' using multiple signals.

    Uses a scoring system combining:
      - Sequence length (chromosomes are typically >500 kb)
      - KNN distance to plasmid training data (chromosomes are distant)
      - Plasmid-typical gene markers in FASTA headers
      - 4-mer profile similarity to known plasmid groups

    Parameters:
        records: list of dicts with 'sequence', 'plasmid_id', 'length', etc.
        vectors: precomputed 4-mer vectors (n_records, 256); computed if None.

    Returns:
        list of dicts with keys:
            plasmid_id, length_bp, classification ('plasmid'/'chromosome'),
            confidence (0-100), reason (explanation string),
            score (raw score, higher = more likely plasmid)
    """
    if not records:
        return []

    # Load training data for distance-based classification
    try:
        ref_data = np.load(CLASSIFIER_PATH, allow_pickle=True)
        X_train = ref_data["X"].astype(np.float64)
        centroids = ref_data["centroids"]
        group_names = list(ref_data["group_names"])
    except Exception:
        X_train = None
        centroids = None
        group_names = []

    # Compute 4-mer vectors if not provided
    if vectors is None:
        seqs = [r["sequence"] for r in records]
        vectors = compute_kmer_vectors(tuple(seqs), k=4)

    results = []
    for i, rec in enumerate(records):
        seq_len = rec["length"]
        plasmid_id = rec["plasmid_id"]
        header = rec.get("plasmid_id", "").lower()
        score = 0      # positive = plasmid, negative = chromosome
        reasons = []

        # ── Signal 1: Sequence length ────────────────────────────────────
        if seq_len > 1_000_000:
            score -= 50
            reasons.append(f"Very large ({seq_len/1e6:.1f} Mb)")
        elif seq_len > CHROMOSOMAL_THRESHOLD:
            score -= 30
            reasons.append(f"Large ({seq_len/1000:.0f} kb, >{CHROMOSOMAL_THRESHOLD/1000:.0f} kb)")
        elif seq_len > 350_000:
            score -= 10
            reasons.append(f"Borderline size ({seq_len/1000:.0f} kb)")
        elif seq_len < 300_000:
            score += 20
            reasons.append(f"Plasmid-range size ({seq_len/1000:.0f} kb)")
        if seq_len < 20_000:
            score += 10
            reasons.append("Small sequence (likely plasmid)")

        # ── Signal 2: Distance to nearest plasmid training vector ────────
        if X_train is not None:
            from scipy.spatial.distance import cdist
            query_vec = vectors[i:i+1]
            dists = cdist(query_vec, X_train, metric="cosine")[0]
            nn_dist = float(np.min(dists))
            mean_dist = float(np.mean(dists))

            # Also compute distance to nearest centroid
            centroid_dists = cdist(query_vec, centroids, metric="cosine")[0]
            min_centroid_dist = float(np.min(centroid_dists))

            if nn_dist < 0.02:
                score += 30
                reasons.append(f"Very close to known plasmid (d={nn_dist:.4f})")
            elif nn_dist < 0.05:
                score += 20
                reasons.append(f"Close to known plasmid (d={nn_dist:.4f})")
            elif nn_dist < 0.10:
                score += 5
                reasons.append(f"Moderate distance to plasmid DB (d={nn_dist:.4f})")
            elif nn_dist < 0.20:
                score -= 10
                reasons.append(f"Distant from plasmid DB (d={nn_dist:.4f})")
            else:
                score -= 25
                reasons.append(f"Very distant from all plasmids (d={nn_dist:.4f})")

            # Centroid distance — chromosomes tend to be far from all centroids
            if min_centroid_dist > 0.25:
                score -= 15
                reasons.append(f"Far from all Inc centroids (d={min_centroid_dist:.4f})")

        # ── Signal 3: FASTA header keywords ──────────────────────────────
        plasmid_keywords = ["plasmid", "unnamed", "p0", "plas"]
        chromo_keywords = ["chromosome", "chr ", "chr1", "chr2", "complete genome",
                           "whole genome", "genomic"]
        if any(kw in header for kw in plasmid_keywords):
            score += 15
            reasons.append("Header suggests plasmid")
        if any(kw in header for kw in chromo_keywords):
            score -= 20
            reasons.append("Header suggests chromosome")

        # ── Signal 4: Inc group confidence ───────────────────────────────
        inc_conf = rec.get("inc_confidence", 0)
        if inc_conf and inc_conf > 60:
            score += 15
            reasons.append(f"Strong Inc group match ({inc_conf:.0f}%)")
        elif inc_conf and inc_conf < 30:
            score -= 5
            reasons.append(f"Weak Inc group match ({inc_conf:.0f}%)")

        # ── Final classification ─────────────────────────────────────────
        # Three categories:
        #   1. "plasmid" (confidence >= 95%) — gets pLIN code + all modules
        #   2. "incomplete_plasmid" (plasmid-like but <95% conf) — no pLIN,
        #      but AMR/mobility/other modules still run
        #   3. "chromosome" — excluded entirely
        confidence = min(99, 50 + abs(score))

        if score <= -10:
            classification = "chromosome"
        elif score >= 10 and confidence >= 95:
            classification = "plasmid"
        else:
            # Borderline or moderate plasmid signal — flag as incomplete
            classification = "incomplete_plasmid"

        results.append({
            "plasmid_id": plasmid_id,
            "length_bp": seq_len,
            "classification": classification,
            "confidence": round(confidence, 1),
            "score": score,
            "reason": "; ".join(reasons),
        })

    return results


def merge_multicontig_plasmids(records):
    """Merge plasmid contigs from the same source file + Inc type into single records.

    When a multi-FASTA assembly file contains multiple plasmid contigs with the
    same Inc type, they are likely fragments of the same replicon. This function
    concatenates them (with 100-N spacer) so they receive a single pLIN code.

    Contigs with different Inc types within the same file are kept separate
    (they represent genuinely different plasmids).

    Returns:
        merged_records: list of dicts (same schema as input records)
        merge_map: dict mapping merged plasmid_id → list of original contig plasmid_ids
    """
    from collections import defaultdict

    # Group records by (source_file, inc_type)
    groups = defaultdict(list)
    for rec in records:
        key = (rec.get("source_file", "unknown"), rec.get("inc_type", "Unknown"))
        groups[key].append(rec)

    merged_records = []
    merge_map = {}

    for (src_file, inc_type), group_recs in groups.items():
        if len(group_recs) == 1:
            # Single contig — pass through unchanged
            rec = group_recs[0]
            merged_records.append(rec)
            merge_map[rec["plasmid_id"]] = [rec["plasmid_id"]]
        else:
            # Multiple contigs from same file with same Inc type — merge
            # Sort by contig index to maintain original order
            group_recs.sort(key=lambda r: r.get("_contig_index", 0))

            # Concatenate sequences with N spacer
            spacer = "N" * 100
            merged_seq = spacer.join(r["sequence"] for r in group_recs)
            original_ids = [r["plasmid_id"] for r in group_recs]

            # Use first contig as template for metadata
            template = group_recs[0].copy()
            merged_id = f"{src_file}|{inc_type}|merged({len(group_recs)})"
            template["plasmid_id"] = merged_id
            template["sequence"] = merged_seq
            template["length"] = len(merged_seq)
            template["_merged_from"] = original_ids
            template["_merged_count"] = len(group_recs)

            # Use best confidence among merged contigs
            confidences = [r.get("inc_confidence", 0) or 0 for r in group_recs]
            if confidences:
                best_idx = confidences.index(max(confidences))
                best_rec = group_recs[best_idx]
                template["inc_confidence"] = best_rec.get("inc_confidence")
                template["inc_best_match"] = best_rec.get("inc_best_match")
                template["inc_top5_candidates"] = best_rec.get("inc_top5_candidates")
                template["inc_is_low_confidence"] = best_rec.get("inc_is_low_confidence")

            # Only skip pLIN if ALL contigs in the group were marked _skip_plin
            if all(r.get("_skip_plin", False) for r in group_recs):
                template["_skip_plin"] = True
            else:
                template.pop("_skip_plin", None)

            merged_records.append(template)
            merge_map[merged_id] = original_ids

    return merged_records, merge_map


# ══════════════════════════════════════════════════════════════════════════════
#  L3: ASSEMBLY COMPLETENESS ASSESSMENT
# ══════════════════════════════════════════════════════════════════════════════

def assess_assembly_completeness(records, prodigal_summary_df=None):
    """Assess assembly completeness for each plasmid.

    Computes metrics: contig count, N50 ratio, circular topology signal,
    coding density (if Prodigal available), and composite completeness score.

    Returns DataFrame with completeness metrics per plasmid.
    """
    results = []
    for rec in records:
        seq = str(rec["sequence"]).upper()
        total_len = len(seq)
        pid = rec["plasmid_id"]

        if total_len == 0:
            results.append({"plasmid_id": pid, "total_length": 0,
                            "n_contigs": 0, "n50_ratio": 0, "circular_signal": False,
                            "coding_density_pct": 0, "n_gaps": 0,
                            "completeness_score": 0, "completeness_status": "POOR"})
            continue

        # 1. Contig count — check for N-gaps (>=10 consecutive Ns = contig break)
        import re
        contigs = re.split(r'N{10,}', seq)
        contigs = [c for c in contigs if len(c) > 0]
        n_contigs = max(len(contigs), 1)
        n_gaps = max(n_contigs - 1, 0)

        # 2. N50 ratio
        contig_lengths = sorted([len(c) for c in contigs], reverse=True)
        cumsum = 0
        n50 = contig_lengths[0]
        half = total_len / 2
        for cl in contig_lengths:
            cumsum += cl
            if cumsum >= half:
                n50 = cl
                break
        n50_ratio = n50 / total_len if total_len > 0 else 0

        # 3. Circular topology signal — check overlap between first and last 500bp
        circular_signal = False
        check_len = min(500, total_len // 4)
        if check_len >= 50:
            head = seq[:check_len]
            tail = seq[-check_len:]
            # Simple identity: count matching bases
            matches = sum(1 for a, b in zip(head, tail) if a == b)
            identity = matches / check_len
            circular_signal = identity >= 0.90

        # 4. Coding density (from Prodigal if available)
        coding_density = 0.0
        if prodigal_summary_df is not None and len(prodigal_summary_df) > 0:
            match = prodigal_summary_df[
                prodigal_summary_df["source_file"].str.contains(pid, na=False)
            ]
            if len(match) > 0:
                coding_density = float(match.iloc[0].get("coding_density_pct", 0))

        # 5. Composite completeness score (0-100)
        score = 0
        score += 40 if n_contigs == 1 else max(0, 40 - (n_contigs - 1) * 10)
        score += 20 if n50_ratio > 0.9 else int(20 * n50_ratio)
        score += 20 if circular_signal else 0
        if coding_density > 0:
            score += 10 if coding_density > 80 else int(10 * coding_density / 80)
        else:
            score += 5  # neutral if no Prodigal data
        score += 10 if n_gaps == 0 else max(0, 10 - n_gaps * 2)
        score = min(100, max(0, score))

        # Status
        if score >= 80:
            status = "COMPLETE"
        elif score >= 60:
            status = "NEAR-COMPLETE"
        elif score >= 40:
            status = "FRAGMENTED"
        else:
            status = "POOR"

        results.append({
            "plasmid_id": pid,
            "total_length": total_len,
            "n_contigs": n_contigs,
            "n50_ratio": round(n50_ratio, 3),
            "circular_signal": circular_signal,
            "coding_density_pct": round(coding_density, 1),
            "n_gaps": n_gaps,
            "completeness_score": score,
            "completeness_status": status,
        })

    return pd.DataFrame(results) if results else pd.DataFrame()


def detect_duplicates(records):
    """Detect identical or near-identical sequences among uploaded files.

    Returns list of (plasmid_a, plasmid_b, cosine_distance) tuples.
    """
    if len(records) < 2:
        return []
    from scipy.spatial.distance import pdist, squareform
    vecs = np.array([_kmer_vector_single(str(r["sequence"])) for r in records])
    dists = squareform(pdist(vecs, metric="cosine"))
    duplicates = []
    seen = set()
    for i in range(len(records)):
        for j in range(i + 1, len(records)):
            if dists[i, j] < 0.0001:
                key = (records[i]["plasmid_id"], records[j]["plasmid_id"])
                if key not in seen:
                    duplicates.append((records[i]["plasmid_id"],
                                       records[j]["plasmid_id"],
                                       round(float(dists[i, j]), 8)))
                    seen.add(key)
    return duplicates


# PlasmidFinder FSA filename → Inc group name mapping
# Covers the FSA files present in ~/plasmidfinder_db/
_PLASMIDFINDER_FSA_TO_INC = {
    "Rep1.fsa":       "IncI1",
    "Rep2.fsa":       "IncI2",
    "Rep3.fsa":       "IncF",
    "Rep7.fsa":       "IncX1",
    "Rep9.fsa":       "IncHI2",
    "Rep10.fsa":      "IncR",
    "Inc18.fsa":      "IncFII",
    "NT_Rep.fsa":     "IncN",
    "Rep_trans.fsa":  "IncHI1",
    "RepL.fsa":       "ColE",
    "Rep3_1.fsa":     "IncFIC",
}

def _replicon_blast_inc_types(sequence: str) -> list:
    """Run BLASTn of sequence against PlasmidFinder replicon database.

    Uses enterobacteriales.fsa (combined FSA with all Inc groups in headers).
    Coverage is measured against the replicon subject sequence (not the query plasmid).
    Returns list of (inc_group, pct_identity, coverage) tuples for hits passing thresholds.
    Only called for sequences > REPLICON_BLAST_SIZE_THRESHOLD bp.
    Returns [] if BLAST unavailable or DB missing.
    """
    import subprocess, tempfile, re

    # Prefer Homebrew blastn, fall back to PATH
    for candidate in ["/opt/homebrew/bin/blastn", "blastn"]:
        try:
            subprocess.run([candidate, "-version"], capture_output=True, timeout=5)
            blastn = candidate
            break
        except Exception:
            blastn = None

    db_dir = PLASMIDFINDER_DB_DIR
    combined_fsa = os.path.join(db_dir, "enterobacteriales.fsa")
    if not blastn or not os.path.isfile(combined_fsa):
        return []

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as qf:
        qf.write(f">query\n{sequence}\n")
        query_path = qf.name

    detected = {}
    try:
        result = subprocess.run(
            [blastn, "-query", query_path, "-subject", combined_fsa,
             "-outfmt", "6 qseqid sseqid pident length slen",
             "-perc_identity", str(REPLICON_BLAST_IDENTITY),
             "-task", "blastn", "-evalue", "1e-5", "-max_hsps", "50"],
            capture_output=True, text=True, timeout=60
        )
        for line in result.stdout.strip().split("\n"):
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            sseqid  = parts[1]   # e.g. "IncHI2_1__BX664015" or "IncHI2A_1__BX664015"
            pident  = float(parts[2])
            aln_len = int(parts[3])
            slen    = int(parts[4])
            # Coverage of the replicon sequence (subject)
            coverage = (aln_len / slen) * 100 if slen > 0 else 0
            if pident < REPLICON_BLAST_IDENTITY or coverage < REPLICON_BLAST_COVERAGE:
                continue
            # Parse Inc group from subject header: "IncHI2A_1__BX664015" → "IncHI2A" → normalise to "IncHI2"
            m = re.match(r"(Inc[A-Za-z0-9]+|Col[A-Za-z0-9]+|rep[A-Za-z0-9]+)", sseqid)
            if not m:
                continue
            raw_group = m.group(1)
            # Collapse sub-types: IncHI2A → IncHI2, IncHI1B → IncHI1, IncN2/IncN3 → IncN
            inc_group = re.sub(r"([A-Z])$", "", raw_group)   # strip trailing letter sub-type
            inc_group = re.sub(r"\d+$", "", inc_group) if inc_group.endswith(("2","3")) and "HI" not in inc_group else inc_group
            # Keep only groups in our 28-group classifier
            if inc_group not in ("IncA","IncAC2","IncC","IncF","IncFIB","IncFIBK","IncFIC","IncFII",
                                  "IncHI1","IncHI2","IncI","IncI1","IncI2","IncN","IncR",
                                  "IncX1","IncX3","IncX4","ColE","ColRNAI"):
                continue
            if inc_group not in detected or detected[inc_group][0] < pident:
                detected[inc_group] = (pident, coverage)
    except Exception:
        pass
    finally:
        try:
            os.unlink(query_path)
        except Exception:
            pass

    return [(grp, pid, cov) for grp, (pid, cov) in sorted(detected.items(), key=lambda x: -x[1][0])]


def classify_inc_group(sequence, group_names, classifier):
    """Classify a plasmid to its Inc group using KNN or centroid distance.

    Returns dict with:
        - predicted_group: str (best Inc type, "Multiple Inc", or "Unknown/Novel")
        - confidence: float (0-1)
        - is_low_confidence: bool
        - is_multiple_inc: bool (True if 2+ Inc types detected above threshold)
        - multiple_inc_types: list of (inc_type, confidence) for detected Inc types
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

        # Check for multiple Inc types (multi-replicon plasmids)
        # Find all Inc types above the multi-Inc threshold
        high_conf_incs = [(inc, conf) for inc, conf in sorted_candidates if conf >= MULTI_INC_THRESHOLD]
        is_multiple_inc = len(high_conf_incs) >= 2

        # Secondary signal: replicon BLAST for large plasmids where KNN gives 100% to one group.
        # KNN k=5 cannot detect a secondary replicon when all 5 neighbours are the same group.
        blast_inc_types = []
        blast_used = False
        if (not is_multiple_inc
                and confidence >= 0.95
                and len(sequence) >= REPLICON_BLAST_SIZE_THRESHOLD
                and os.path.isdir(PLASMIDFINDER_DB_DIR)):
            blast_hits = _replicon_blast_inc_types(sequence)
            # Only promote to multi-replicon if BLAST finds a second Inc group
            # that is different from the KNN winner and passes thresholds
            blast_others = [g for g, pid, cov in blast_hits if g != best_group]
            if blast_others:
                blast_inc_types = blast_hits
                blast_used = True
                # Merge: KNN winner + BLAST additional groups
                blast_groups = [g for g, _, _ in blast_hits]
                combined = [best_group] + [g for g in blast_groups if g != best_group]
                high_conf_incs = [(g, confidence if g == best_group else 0.0) for g in combined]
                is_multiple_inc = True

        # Determine the predicted group label
        if is_low_confidence:
            predicted_group = "Unknown/Novel"
        elif is_multiple_inc:
            inc_names = [inc for inc, _ in high_conf_incs[:3]]
            predicted_group = "Multiple: " + "/".join(inc_names)
        else:
            predicted_group = best_group

        return {
            "predicted_group": predicted_group,
            "best_match": best_group,
            "confidence": confidence,
            "is_low_confidence": is_low_confidence,
            "is_multiple_inc": is_multiple_inc,
            "multiple_inc_types": high_conf_incs if is_multiple_inc else [],
            "top5_candidates": sorted_candidates,
            "all_probabilities": proba_dict,
            "blast_replicon_hits": blast_inc_types,
            "blast_used": blast_used,
            # CV-based quality metadata (populated by caller if available)
            "cv_f1": None,
            "is_confusion_pair": False,
            "confusion_note": "",
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

        # Check for multiple Inc types
        high_conf_incs = [(inc, conf) for inc, conf in sorted_candidates if conf >= MULTI_INC_THRESHOLD]
        is_multiple_inc = len(high_conf_incs) >= 2

        # Determine the predicted group label
        if is_low_confidence:
            predicted_group = "Unknown/Novel"
        elif is_multiple_inc:
            inc_names = [inc for inc, _ in high_conf_incs[:3]]
            predicted_group = "Multiple: " + "/".join(inc_names)
        else:
            predicted_group = best_group

        return {
            "predicted_group": predicted_group,
            "best_match": best_group,
            "confidence": confidence,
            "is_low_confidence": is_low_confidence,
            "is_multiple_inc": is_multiple_inc,
            "multiple_inc_types": high_conf_incs if is_multiple_inc else [],
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
    to the hierarchical levels (L1→L6).

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
        "A": 0.99,   # L1 — nearly all within-group distances below this
        "B": 0.95,   # L2
        "C": 0.75,   # L3
        "D": 0.50,   # L4 — median distance
        "E": 0.25,   # L5
        "F": 0.05,   # L6 — only very close pairs
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
#  L4: DATABASE COVERAGE & NOVELTY DETECTION
# ══════════════════════════════════════════════════════════════════════════════

def assess_database_coverage(query_vectors, query_inc_types, query_nn_distances=None):
    """Assess how well the reference database covers query plasmids.

    Computes: NN distance percentile within Inc group, novelty flags,
    database representation scores.

    Returns DataFrame with coverage metrics per plasmid.
    """
    if not os.path.exists(CLASSIFIER_PATH):
        return pd.DataFrame()

    data = np.load(CLASSIFIER_PATH, allow_pickle=True)
    X_train = data["X"]
    y_train = data["y"]
    group_names = [str(g) for g in data["group_names"]]

    # Pre-compute within-group distance distributions
    from scipy.spatial.distance import cdist
    group_stats = {}
    for i, name in enumerate(group_names):
        mask = y_train == i
        n_samples = int(mask.sum())
        group_stats[name] = {"n_training": n_samples}

        if n_samples >= 10:
            X_group = X_train[mask]
            # Sample for efficiency
            if n_samples > 200:
                rng = np.random.default_rng(42)
                idx = rng.choice(n_samples, 200, replace=False)
                X_sub = X_group[idx]
            else:
                X_sub = X_group
            dists = pdist(X_sub.astype(np.float64), metric="cosine")
            group_stats[name]["p50"] = float(np.percentile(dists, 50))
            group_stats[name]["p95"] = float(np.percentile(dists, 95))
            group_stats[name]["p99"] = float(np.percentile(dists, 99))

    # Compute centroid distances for each query
    centroids = {}
    for i, name in enumerate(group_names):
        mask = y_train == i
        if mask.sum() > 0:
            centroids[name] = X_train[mask].mean(axis=0)

    results = []
    for j in range(len(query_inc_types)):
        inc = query_inc_types[j] if j < len(query_inc_types) else "Unknown"
        nn_dist = query_nn_distances[j] if query_nn_distances is not None and j < len(query_nn_distances) else None

        stats = group_stats.get(inc, {})
        n_training = stats.get("n_training", 0)

        # Representation flag
        if n_training >= 50:
            representation = "Well-represented"
        elif n_training >= 20:
            representation = "Moderate"
        elif n_training > 0:
            representation = "Under-represented"
        else:
            representation = "Not in database"

        # Novelty detection
        novelty_flag = "None"
        nn_percentile = None
        centroid_dist = None

        if inc in centroids and j < len(query_vectors):
            cd = float(1 - np.dot(query_vectors[j], centroids[inc]) /
                       (np.linalg.norm(query_vectors[j]) * np.linalg.norm(centroids[inc]) + 1e-10))
            centroid_dist = round(cd, 6)

        if nn_dist is not None:
            if nn_dist > 0.050:
                novelty_flag = "Potentially novel lineage"
            elif nn_dist > 0.020:
                novelty_flag = "Divergent"

            # Percentile within group
            p99 = stats.get("p99")
            p95 = stats.get("p95")
            p50 = stats.get("p50")
            if p99 is not None:
                if nn_dist > p99:
                    nn_percentile = 99
                elif nn_dist > p95:
                    nn_percentile = 95
                elif nn_dist > p50:
                    nn_percentile = 75
                else:
                    nn_percentile = 50

        # Traffic light
        if novelty_flag == "Potentially novel lineage" or representation == "Not in database":
            coverage_indicator = "RED"
        elif novelty_flag == "Divergent" or representation == "Under-represented":
            coverage_indicator = "YELLOW"
        else:
            coverage_indicator = "GREEN"

        results.append({
            "inc_type": inc,
            "n_training_samples": n_training,
            "representation": representation,
            "nn_distance": round(nn_dist, 6) if nn_dist is not None else None,
            "nn_percentile": nn_percentile,
            "centroid_distance": centroid_dist,
            "novelty_flag": novelty_flag,
            "coverage_indicator": coverage_indicator,
        })

    return pd.DataFrame(results) if results else pd.DataFrame()


# ══════════════════════════════════════════════════════════════════════════════
#  L6: NOVEL INC GROUP DISCOVERY
# ══════════════════════════════════════════════════════════════════════════════

def discover_novel_inc_groups(vectors, inc_predictions, confidence_scores, records):
    """Cluster plasmids with Unknown/Novel Inc classification to discover
    putative new Inc/rep type groups.

    Returns DataFrame with novel group assignments and nearest known Inc group.
    """
    if not os.path.exists(CLASSIFIER_PATH):
        return pd.DataFrame()

    # Find novel plasmids (confidence < 40%)
    novel_indices = [i for i, conf in enumerate(confidence_scores) if conf < 40]
    if len(novel_indices) < 3:
        return pd.DataFrame()

    novel_vectors = vectors[novel_indices]
    novel_ids = [records[i]["plasmid_id"] for i in novel_indices]

    # Cluster novel plasmids
    from scipy.cluster.hierarchy import linkage, fcluster
    dist_matrix = pdist(novel_vectors.astype(np.float64), metric="cosine")
    Z = linkage(dist_matrix, method="average")
    clusters = fcluster(Z, t=0.050, criterion="distance")  # L3 threshold

    # Load centroids for distance comparison
    data = np.load(CLASSIFIER_PATH, allow_pickle=True)
    group_names = [str(g) for g in data["group_names"]]
    X_train = data["X"]
    y_train = data["y"]

    centroids = {}
    for i, name in enumerate(group_names):
        mask = y_train == i
        if mask.sum() > 0:
            centroids[name] = X_train[mask].mean(axis=0)

    # Analyze each cluster
    from collections import Counter
    cluster_counts = Counter(clusters)
    results = []

    for idx, (novel_idx, cluster_id) in enumerate(zip(novel_indices, clusters)):
        n_members = cluster_counts[cluster_id]
        pid = novel_ids[idx]

        # Find nearest known Inc group
        nearest_inc = "Unknown"
        nearest_dist = 1.0
        query_vec = novel_vectors[idx]
        for inc_name, centroid in centroids.items():
            d = float(1 - np.dot(query_vec, centroid) /
                      (np.linalg.norm(query_vec) * np.linalg.norm(centroid) + 1e-10))
            if d < nearest_dist:
                nearest_dist = d
                nearest_inc = inc_name

        is_putative_group = n_members >= 3

        results.append({
            "plasmid_id": pid,
            "novel_cluster": int(cluster_id),
            "cluster_size": n_members,
            "is_putative_novel_group": is_putative_group,
            "nearest_known_inc": nearest_inc,
            "distance_to_nearest": round(nearest_dist, 4),
            "confidence": round(confidence_scores[novel_indices[idx]], 1),
        })

    return pd.DataFrame(results) if results else pd.DataFrame()


# ══════════════════════════════════════════════════════════════════════════════
#  MOBILITY PREDICTION
# ══════════════════════════════════════════════════════════════════════════════

def classify_mobility(amr_df, source_file, mobsuite_df=None):
    """Classify plasmid mobility using best available data.

    Priority: MOBsuite (gold standard) > AMRFinderPlus gene scan > Non-mobilizable.

    Categories:
    - Conjugative: has tra/trb/virB/trw transfer genes or MOBsuite says conjugative
    - Mobilizable: has mob/nik genes or MOBsuite says mobilizable (needs helper)
    - Non-mobilizable: no detectable transfer/mobilization genes

    Returns dict with keys:
        mobility, genes, detail, source, relaxase_family, mpf_type
    """
    sf = source_file.replace(".fasta", "").replace(".fa", "").replace(".fna", "")

    # ── Tier 1: MOBsuite (if available) ──
    if mobsuite_df is not None and len(mobsuite_df) > 0:
        mob_hits = mobsuite_df[mobsuite_df["source_file"] == sf]
        if len(mob_hits) > 0:
            row = mob_hits.iloc[0]
            mobility = str(row.get("predicted_mobility", "")).lower()
            relaxase = str(row.get("relaxase_type(s)", ""))
            mpf = str(row.get("mpf_type", ""))
            if relaxase in ("-", "nan", ""):
                relaxase = ""
            if mpf in ("-", "nan", ""):
                mpf = ""
            detail_parts = []
            if relaxase:
                detail_parts.append(f"Relaxase: {relaxase}")
            if mpf:
                detail_parts.append(f"MPF: {mpf}")
            detail = "; ".join(detail_parts) if detail_parts else "MOBsuite classification"

            if "conjugative" in mobility:
                return {
                    "mobility": "Conjugative",
                    "genes": [g for g in [relaxase, mpf] if g],
                    "detail": detail, "source": "MOBsuite",
                    "relaxase_family": relaxase, "mpf_type": mpf,
                }
            elif "mobilizable" in mobility:
                return {
                    "mobility": "Mobilizable",
                    "genes": [relaxase] if relaxase else [],
                    "detail": detail, "source": "MOBsuite",
                    "relaxase_family": relaxase, "mpf_type": mpf,
                }
            else:
                return {
                    "mobility": "Non-mobilizable",
                    "genes": [], "detail": detail, "source": "MOBsuite",
                    "relaxase_family": "", "mpf_type": "",
                }

    # ── Tier 2: AMRFinderPlus gene scan (expanded) ──
    if amr_df is None or len(amr_df) == 0:
        return {
            "mobility": "Unknown", "genes": [], "detail": "No AMR data",
            "source": "None", "relaxase_family": "", "mpf_type": "",
        }

    hits = amr_df[amr_df["source_file"] == sf]
    if len(hits) == 0:
        return {
            "mobility": "Non-mobilizable", "genes": [],
            "detail": "No genes detected", "source": "AMRFinderPlus",
            "relaxase_family": "", "mpf_type": "",
        }

    gene_col = "Element symbol" if "Element symbol" in hits.columns else None
    name_col = "Element name" if "Element name" in hits.columns else None

    if gene_col is None:
        return {
            "mobility": "Unknown", "genes": [],
            "detail": "No gene symbol column", "source": "AMRFinderPlus",
            "relaxase_family": "", "mpf_type": "",
        }

    all_genes = hits[gene_col].tolist()
    all_names = hits[name_col].tolist() if name_col else []

    # Check for conjugative markers (prefix match)
    conj_found = []
    for gene in all_genes:
        gene_str = str(gene)
        for prefix in MOBILITY_PREFIXES_CONJUGATIVE:
            if gene_str.startswith(prefix):
                conj_found.append(gene_str)
                break

    # Scan gene names/descriptions for conjugative keywords
    for i, name in enumerate(all_names):
        name_lower = str(name).lower()
        for keyword in MOBILITY_KEYWORDS_CONJUGATIVE:
            if keyword in name_lower:
                g = str(all_genes[i])
                if g not in conj_found:
                    conj_found.append(g)
                break

    # Check for mobilizable markers (prefix match)
    mob_found = []
    for gene in all_genes:
        gene_str = str(gene)
        for prefix in MOBILITY_PREFIXES_MOBILIZABLE:
            if gene_str.startswith(prefix):
                mob_found.append(gene_str)
                break

    # Scan gene names/descriptions for mobilizable keywords
    for i, name in enumerate(all_names):
        name_lower = str(name).lower()
        for keyword in MOBILITY_KEYWORDS_MOBILIZABLE:
            if keyword in name_lower:
                g = str(all_genes[i])
                if g not in mob_found and g not in conj_found:
                    mob_found.append(g)
                break

    if conj_found:
        unique_genes = list(set(conj_found))
        return {
            "mobility": "Conjugative", "genes": unique_genes,
            "detail": f"{len(unique_genes)} transfer gene(s)",
            "source": "AMRFinderPlus", "relaxase_family": "", "mpf_type": "",
        }
    elif mob_found:
        unique_genes = list(set(mob_found))
        return {
            "mobility": "Mobilizable", "genes": unique_genes,
            "detail": f"{len(unique_genes)} mobilization gene(s)",
            "source": "AMRFinderPlus", "relaxase_family": "", "mpf_type": "",
        }
    else:
        return {
            "mobility": "Non-mobilizable", "genes": [],
            "detail": "No transfer/mobilization genes detected",
            "source": "AMRFinderPlus", "relaxase_family": "", "mpf_type": "",
        }


# ══════════════════════════════════════════════════════════════════════════════
#  L10: MOBILE GENETIC ELEMENT BOUNDARY DETECTION
# ══════════════════════════════════════════════════════════════════════════════

# MGE marker gene patterns
_MGE_PATTERNS = {
    "IS_element": [
        r'\bIS\d+', r'\bIS[A-Z][a-z]+\d*', r'\btnp[A-Z]?\b', r'\btransposase\b',
        r'\bIS\b', r'\binsertion.sequence\b',
    ],
    "integrase": [
        r'\bintegrase\b', r'\bintI\b', r'\bphage.integrase\b', r'\bsite.specific.recombinase\b',
    ],
    "recombinase": [
        r'\brecombinase\b', r'\bresolvase\b', r'\binvertase\b',
    ],
    "excisionase": [
        r'\bexcisionase\b', r'\bxis\b',
    ],
}

_BACKBONE_PATTERNS = [
    r'\brep[A-Z]?\b', r'\breplication\b', r'\bpar[A-Z]?\b', r'\bpartition\b',
    r'\btoxin\b', r'\bantitoxin\b', r'\bstb[A-Z]?\b', r'\brelB\b', r'\brelE\b',
    r'\bori[TV]?\b',
]

_RESISTANCE_PATTERNS = [
    r'\bbla[A-Z]', r'\baac\b', r'\baph\b', r'\bant\b', r'\berm\b', r'\btet\b',
    r'\bsul\b', r'\bdfr\b', r'\bmcr\b', r'\bqnr\b', r'\bcfr\b', r'\bfos\b',
    r'\bvan[A-Z]\b', r'\bmef\b', r'\bmph\b',
]

_MOBILITY_PATTERNS = [
    r'\btra[A-Z]\b', r'\btrb[A-Z]\b', r'\bmob[A-Z]\b', r'\bvirB\b',
    r'\bnik[A-Z]\b', r'\brelaxase\b', r'\bconjugal\b', r'\bT4SS\b',
]


def detect_mge_boundaries(genes_df, amr_genes=None):
    """Detect mobile genetic element boundaries from gene annotations.

    Scans Prodigal gene names for IS elements, integrases, recombinases,
    and identifies putative composite transposons (IS pairs flanking cargo).

    Returns DataFrame with per-gene MGE classification and region annotations.
    """
    import re

    if genes_df is None or len(genes_df) == 0:
        return pd.DataFrame()

    results = []
    for _, gene in genes_df.iterrows():
        gene_name = str(gene.get("gene_name", "")).lower()
        gene_product = str(gene.get("product", "")).lower()
        combined = f"{gene_name} {gene_product}"

        # Classify gene type
        gene_type = "hypothetical"
        mge_subtype = None

        # Check MGE patterns
        for mge_type, patterns in _MGE_PATTERNS.items():
            for pat in patterns:
                if re.search(pat, combined, re.IGNORECASE):
                    gene_type = "MGE"
                    mge_subtype = mge_type
                    break
            if gene_type == "MGE":
                break

        # Check backbone
        if gene_type == "hypothetical":
            for pat in _BACKBONE_PATTERNS:
                if re.search(pat, combined, re.IGNORECASE):
                    gene_type = "backbone"
                    break

        # Check resistance
        if gene_type == "hypothetical":
            for pat in _RESISTANCE_PATTERNS:
                if re.search(pat, combined, re.IGNORECASE):
                    gene_type = "resistance"
                    break

        # Check mobility
        if gene_type == "hypothetical":
            for pat in _MOBILITY_PATTERNS:
                if re.search(pat, combined, re.IGNORECASE):
                    gene_type = "mobility"
                    break

        results.append({
            "source_file": gene.get("source_file", ""),
            "gene_id": gene.get("gene_id", ""),
            "start": gene.get("start", 0),
            "end": gene.get("end", 0),
            "strand": gene.get("strand", "+"),
            "gene_name": gene.get("gene_name", ""),
            "gene_type": gene_type,
            "mge_subtype": mge_subtype,
        })

    mge_df = pd.DataFrame(results)

    # Detect composite transposons: IS elements flanking cargo genes
    if len(mge_df) > 0:
        for source, group in mge_df.groupby("source_file"):
            is_positions = group[group["mge_subtype"] == "IS_element"].sort_values("start")
            if len(is_positions) >= 2:
                for i in range(len(is_positions) - 1):
                    is1_end = is_positions.iloc[i]["end"]
                    is2_start = is_positions.iloc[i + 1]["start"]
                    # Check if there are cargo genes between IS elements
                    cargo = group[(group["start"] >= is1_end) & (group["end"] <= is2_start)]
                    has_resistance = (cargo["gene_type"] == "resistance").any()
                    if len(cargo) > 0 and has_resistance:
                        # Mark as composite transposon region
                        mask = ((mge_df["source_file"] == source) &
                                (mge_df["start"] >= is_positions.iloc[i]["start"]) &
                                (mge_df["end"] <= is_positions.iloc[i + 1]["end"]))
                        mge_df.loc[mask & (mge_df["gene_type"] == "hypothetical"), "gene_type"] = "MGE-cargo"

    return mge_df


def draw_plasmid_gene_map(mge_df, plasmid_length, plasmid_id):
    """Draw linear gene map with color-coded regions.

    Returns matplotlib figure showing gene architecture with:
    - Blue: backbone genes
    - Red: resistance genes
    - Green: mobility genes
    - Yellow: IS/MGE elements
    - Orange: MGE cargo
    - Gray: hypothetical
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    color_map = {
        "backbone": "#1E88E5",
        "resistance": "#E53935",
        "mobility": "#43A047",
        "MGE": "#FFC107",
        "MGE-cargo": "#FB8C00",
        "hypothetical": "#BDBDBD",
    }

    fig, ax = plt.subplots(1, 1, figsize=(14, 2.5))

    # Draw backbone line
    ax.plot([0, plasmid_length], [0, 0], color="#E0E0E0", linewidth=8, solid_capstyle="round")

    # Draw genes as arrows
    for _, gene in mge_df.iterrows():
        start = gene["start"]
        end = gene["end"]
        strand = gene.get("strand", "+")
        gene_type = gene["gene_type"]
        color = color_map.get(gene_type, "#BDBDBD")

        width = end - start
        direction = 1 if strand == "+" else -1
        height = 0.4

        # Arrow glyph
        if direction > 0:
            arrow = mpatches.FancyArrow(start, 0, width, 0, width=height,
                                         head_width=height * 1.3, head_length=min(width * 0.2, 500),
                                         fc=color, ec="white", linewidth=0.5)
        else:
            arrow = mpatches.FancyArrow(end, 0, -width, 0, width=height,
                                         head_width=height * 1.3, head_length=min(width * 0.2, 500),
                                         fc=color, ec="white", linewidth=0.5)
        ax.add_patch(arrow)

    ax.set_xlim(-plasmid_length * 0.02, plasmid_length * 1.02)
    ax.set_ylim(-1.2, 1.5)
    ax.set_xlabel("Position (bp)", fontsize=9)
    ax.set_title(f"{plasmid_id} — Gene Architecture ({plasmid_length:,} bp)", fontsize=11)
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Legend
    legend_patches = [mpatches.Patch(color=c, label=l.replace("MGE-cargo", "MGE cargo"))
                      for l, c in color_map.items() if l != "MGE-cargo"]
    legend_patches.append(mpatches.Patch(color=color_map["MGE-cargo"], label="MGE cargo"))
    ax.legend(handles=legend_patches, loc="upper right", fontsize=7, ncol=3,
              framealpha=0.8)

    plt.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  WITHIN-GROUP SEQUENCE ALIGNMENT
# ══════════════════════════════════════════════════════════════════════════════
#
# A pLIN code (or a shared code prefix) says two plasmids are compositionally
# similar, but composition similarity is not itself evidence of which regions
# are actually shared. This module runs pairwise blastn between members of a
# user-selected pLIN group and reports coverage, identity, and the coordinates
# of shared vs. non-shared sequence directly, plus which AMR genes (if any)
# fall inside vs. outside the shared blocks — the same analysis performed by
# hand for Supplementary Table S5 during peer review, now a reusable feature
# rather than a one-off validation script.

def run_within_group_alignment(records, plasmid_ids, blastn_binary=None, min_align_len=500):
    """Pairwise blastn alignment between all selected plasmids.

    Parameters
    ----------
    records : list of dict
        Plasmid records as stored in st.session_state.records (must include
        "plasmid_id" and "sequence").
    plasmid_ids : list of str
        Subset of plasmid_ids to align, all pairwise combinations.
    blastn_binary : str, optional
        Path to blastn; auto-detected via detect_blastn() if not given.
    min_align_len : int
        Discard HSPs shorter than this (bp) as noise.

    Returns
    -------
    dict with keys:
        "pairs": DataFrame, one row per plasmid pair — query/subject id,
            length of each, total length covered, % coverage (of the
            shorter sequence), weighted % identity, n_blocks.
        "blocks": DataFrame, one row per aligned block (HSP) — pair,
            query/subject coordinates, length, % identity. This is the
            data needed to say *which* regions are shared, not just how
            much.
        "error": str or None. If blastn/makeblastdb are unavailable, or a
            pair fails, this is set and "pairs"/"blocks" may be empty.
    """
    if blastn_binary is None:
        blastn_binary, makeblastdb_binary = detect_blastn()
    else:
        _, makeblastdb_binary = detect_blastn()

    if not blastn_binary or not makeblastdb_binary:
        return {"pairs": pd.DataFrame(), "blocks": pd.DataFrame(),
                "error": "BLAST+ (blastn/makeblastdb) not found. Install BLAST+ to use "
                         "within-group alignment (see Installation docs)."}

    # plasmid_id is not guaranteed unique across records — e.g. multi-contig
    # assemblies with generic contig names ("1", "unnamed", ...) can collide
    # across different uploaded files. Silently keying a dict by plasmid_id
    # would drop one and silently align the wrong sequence, so check first.
    selected_records = [r for r in records if r["plasmid_id"] in plasmid_ids]
    id_counts = {}
    for r in selected_records:
        id_counts[r["plasmid_id"]] = id_counts.get(r["plasmid_id"], 0) + 1
    duplicated = sorted(pid for pid, n in id_counts.items() if n > 1)
    if duplicated:
        return {"pairs": pd.DataFrame(), "blocks": pd.DataFrame(),
                "error": "Duplicate plasmid_id among selected records (cannot tell them "
                         f"apart): {', '.join(duplicated)}. This can happen with generic "
                         "contig names shared across different uploaded files."}

    seq_by_id = {r["plasmid_id"]: r["sequence"] for r in selected_records}
    missing = [pid for pid in plasmid_ids if pid not in seq_by_id]
    if missing:
        return {"pairs": pd.DataFrame(), "blocks": pd.DataFrame(),
                "error": f"Sequence not found in session for: {', '.join(missing)}"}
    if len(seq_by_id) < 2:
        return {"pairs": pd.DataFrame(), "blocks": pd.DataFrame(),
                "error": "Select at least 2 plasmids to align."}

    pair_rows = []
    block_rows = []
    work_dir = tempfile.mkdtemp(prefix="plin_align_")
    try:
        ordered_ids = [pid for pid in plasmid_ids if pid in seq_by_id]
        for i in range(len(ordered_ids)):
            for j in range(i + 1, len(ordered_ids)):
                qid, sid = ordered_ids[i], ordered_ids[j]
                qseq, sseq = seq_by_id[qid], seq_by_id[sid]

                query_path = os.path.join(work_dir, "query.fasta")
                subject_path = os.path.join(work_dir, "subject.fasta")
                with open(query_path, "w") as fh:
                    fh.write(f">{qid}\n{qseq}\n")
                with open(subject_path, "w") as fh:
                    fh.write(f">{sid}\n{sseq}\n")

                try:
                    result = subprocess.run(
                        [blastn_binary, "-query", query_path, "-subject", subject_path,
                         "-outfmt", "6 qstart qend sstart send length pident evalue",
                         "-task", "blastn", "-evalue", "1e-10"],
                        capture_output=True, text=True, timeout=120)
                except Exception as exc:
                    pair_rows.append({
                        "plasmid_1": qid, "plasmid_2": sid,
                        "length_1_bp": len(qseq), "length_2_bp": len(sseq),
                        "coverage_pct": None, "weighted_identity_pct": None,
                        "n_blocks": 0, "error": str(exc),
                    })
                    continue

                hsps = []
                for line in result.stdout.strip().split("\n"):
                    if not line:
                        continue
                    qstart, qend, sstart, send, length, pident, evalue = line.split("\t")
                    length = int(length)
                    if length < min_align_len:
                        continue
                    hsps.append({
                        "plasmid_1": qid, "plasmid_2": sid,
                        "q_start": min(int(qstart), int(qend)), "q_end": max(int(qstart), int(qend)),
                        "s_start": min(int(sstart), int(send)), "s_end": max(int(sstart), int(send)),
                        "length_bp": length, "pct_identity": float(pident),
                    })
                hsps.sort(key=lambda h: h["q_start"])
                block_rows.extend(hsps)

                # Coverage is computed separately on each side (union of that
                # side's aligned intervals, overlaps merged so duplicate/
                # overlapping HSPs aren't double-counted, as a fraction of
                # that side's own length) and then reported as coverage of
                # the shorter sequence specifically — dividing a query-side
                # union by the subject's length (or vice versa) is a
                # mismatched numerator/denominator and can exceed 100% when
                # the query is not the shorter sequence.
                def _merged_union(pairs):
                    if not pairs:
                        return 0
                    pairs = sorted(pairs)
                    merged = [pairs[0]]
                    for s, e in pairs[1:]:
                        if s <= merged[-1][1]:
                            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                        else:
                            merged.append((s, e))
                    return sum(e - s for s, e in merged)

                weighted_identity_num = 0.0
                if hsps:
                    q_covered = _merged_union([(h["q_start"], h["q_end"]) for h in hsps])
                    s_covered = _merged_union([(h["s_start"], h["s_end"]) for h in hsps])
                    covered_of_shorter = q_covered if len(qseq) <= len(sseq) else s_covered
                    weighted_identity_num = sum(h["length_bp"] * h["pct_identity"] for h in hsps)
                else:
                    covered_of_shorter = 0

                shorter_len = min(len(qseq), len(sseq))
                coverage_pct = round(100 * covered_of_shorter / shorter_len, 1) if shorter_len else 0.0
                total_aligned_bp = sum(h["length_bp"] for h in hsps)
                weighted_identity = round(weighted_identity_num / total_aligned_bp, 2) if total_aligned_bp else None

                pair_rows.append({
                    "plasmid_1": qid, "plasmid_2": sid,
                    "length_1_bp": len(qseq), "length_2_bp": len(sseq),
                    "coverage_pct": coverage_pct,
                    "weighted_identity_pct": weighted_identity,
                    "n_blocks": len(hsps), "error": None,
                })
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    return {"pairs": pd.DataFrame(pair_rows), "blocks": pd.DataFrame(block_rows), "error": None}


def annotate_alignment_blocks_with_amr(blocks_df, amr_df):
    """Flag which AMR genes fall inside vs. outside each aligned shared block.

    For every (plasmid_1, plasmid_2) pair in blocks_df, checks amr_df for
    genes on either plasmid and reports whether each gene's coordinates fall
    within any shared block on that plasmid's side of the alignment. This is
    what lets a user see, directly, whether AMR content differs between two
    plasmids that otherwise look nearly identical by pLIN code — the
    question Reviewer 2 asked for by name.

    Returns a DataFrame: plasmid_id, gene, class, start, end, in_shared_block
    (bool), paired_with (the other plasmid_id in the comparison).
    """
    if blocks_df is None or len(blocks_df) == 0 or amr_df is None or len(amr_df) == 0:
        return pd.DataFrame()
    if "Start" not in amr_df.columns or "source_file" not in amr_df.columns:
        return pd.DataFrame()

    rows = []
    for (p1, p2), pair_blocks in blocks_df.groupby(["plasmid_1", "plasmid_2"]):
        q_intervals = [(b["q_start"], b["q_end"]) for _, b in pair_blocks.iterrows()]
        s_intervals = [(b["s_start"], b["s_end"]) for _, b in pair_blocks.iterrows()]

        for plasmid_id, intervals, other in ((p1, q_intervals, p2), (p2, s_intervals, p1)):
            genes = amr_df[amr_df["source_file"] == plasmid_id]
            for _, gene in genes.iterrows():
                try:
                    gstart, gend = int(gene["Start"]), int(gene["Stop"])
                except (ValueError, TypeError, KeyError):
                    continue
                gstart, gend = min(gstart, gend), max(gstart, gend)
                in_block = any(gstart >= s and gend <= e for s, e in intervals) or \
                           any(min(gend, e) - max(gstart, s) > 0 for s, e in intervals)
                rows.append({
                    "plasmid_id": plasmid_id,
                    "paired_with": other,
                    "gene": gene.get("Element symbol", ""),
                    "class": gene.get("Class", ""),
                    "start": gstart, "end": gend,
                    "in_shared_block": bool(in_block),
                })
    return pd.DataFrame(rows)


def draw_alignment_comparison(blocks_row_group, plasmid_1, plasmid_2, len_1, len_2,
                               amr_annotated=None):
    """Draw a two-track linear comparison of one plasmid pair with shared
    blocks connected between tracks, similar in spirit to an Easyfig/ACT
    view but lightweight (matplotlib only, no external alignment viewer).
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.path import Path as MplPath

    fig, ax = plt.subplots(1, 1, figsize=(14, 3.2))
    y_top, y_bot = 1.0, 0.0
    track_h = 0.15

    ax.add_patch(mpatches.Rectangle((0, y_top - track_h / 2), len_1, track_h,
                                     fc="#E0E0E0", ec="none"))
    ax.add_patch(mpatches.Rectangle((0, y_bot - track_h / 2), len_2, track_h,
                                     fc="#E0E0E0", ec="none"))

    cmap = plt.get_cmap("viridis")
    n_blocks = len(blocks_row_group)
    for idx, (_, b) in enumerate(blocks_row_group.iterrows()):
        color = cmap(idx / max(n_blocks - 1, 1)) if n_blocks > 1 else cmap(0.5)
        alpha = min(0.3 + b["pct_identity"] / 200, 0.85)

        ax.add_patch(mpatches.Rectangle((b["q_start"], y_top - track_h / 2),
                                         b["q_end"] - b["q_start"], track_h,
                                         fc=color, ec="none"))
        ax.add_patch(mpatches.Rectangle((b["s_start"], y_bot - track_h / 2),
                                         b["s_end"] - b["s_start"], track_h,
                                         fc=color, ec="none"))

        verts = [(b["q_start"], y_top - track_h / 2), (b["s_start"], y_bot + track_h / 2),
                 (b["s_end"], y_bot + track_h / 2), (b["q_end"], y_top - track_h / 2),
                 (b["q_start"], y_top - track_h / 2)]
        codes = [MplPath.MOVETO, MplPath.LINETO, MplPath.LINETO, MplPath.LINETO, MplPath.CLOSEPOLY]
        ax.add_patch(mpatches.PathPatch(MplPath(verts, codes), fc=color, ec="none", alpha=alpha))

    if amr_annotated is not None and len(amr_annotated) > 0:
        for plasmid_id, y, length in ((plasmid_1, y_top, len_1), (plasmid_2, y_bot, len_2)):
            genes = amr_annotated[amr_annotated["plasmid_id"] == plasmid_id]
            for _, g in genes.iterrows():
                color = "#43A047" if g["in_shared_block"] else "#E53935"
                y_off = track_h if y == y_top else -track_h
                ax.plot([g["start"], g["end"]], [y + y_off, y + y_off],
                        color=color, linewidth=3, solid_capstyle="butt")

    ax.set_xlim(-max(len_1, len_2) * 0.02, max(len_1, len_2) * 1.02)
    ax.set_ylim(-0.5, 1.5)
    ax.set_yticks([y_bot, y_top])
    ax.set_yticklabels([f"{plasmid_2}\n({len_2:,} bp)", f"{plasmid_1}\n({len_1:,} bp)"], fontsize=8)
    ax.set_xlabel("Position (bp)", fontsize=9)
    ax.set_title(f"{plasmid_1} vs {plasmid_2} — shared blocks (color = block; darker = higher identity)",
                 fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if amr_annotated is not None and len(amr_annotated) > 0:
        legend_patches = [
            mpatches.Patch(color="#43A047", label="AMR gene — inside shared block"),
            mpatches.Patch(color="#E53935", label="AMR gene — outside shared block"),
        ]
        ax.legend(handles=legend_patches, loc="upper right", fontsize=7, framealpha=0.8)

    plt.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  OUTBREAK / CLONE DETECTION
# ══════════════════════════════════════════════════════════════════════════════

def detect_outbreak_clusters(plin_df, integrated_df):
    """Flag potential outbreak clusters: plasmids sharing the same full pLIN
    strain-level code (i.e. identical down to L6). Matching AMR resistance
    profile is attached as supplementary evidence when available, but is not
    required — a shared pLIN code alone is sufficient evidence of a likely
    shared plasmid, since pLIN already reflects whole-plasmid similarity.

    Returns list of dicts with cluster info.
    """
    if plin_df is None or len(plin_df) == 0 or "pLIN" not in plin_df.columns:
        return []

    df = plin_df.copy()

    # Merge in AMR gene calls when available, purely for supplementary evidence.
    if integrated_df is not None and len(integrated_df) > 0 and "AMR_genes" in integrated_df.columns:
        amr_lookup = integrated_df.groupby("plasmid_id")["AMR_genes"].apply(
            lambda s: "; ".join(sorted({g.strip() for x in s.fillna("") for g in x.split(";") if g.strip()}))
        )
        df["_amr_genes_str"] = df["plasmid_id"].map(amr_lookup).fillna("")
    else:
        df["_amr_genes_str"] = ""

    clusters = []
    for plin_code, group in df.groupby("pLIN"):
        if len(group) < 2 or not isinstance(plin_code, str):
            continue
        if plin_code == "UNCLASSIFIED" or "NEW" in plin_code:
            continue  # Divergent/novel codes are not a shared-plasmid match

        plasmids = group["plasmid_id"].tolist()
        amr_genes = sorted({
            g.strip() for x in group["_amr_genes_str"] for g in x.split(";") if g.strip()
        })
        amr_consistent = group["_amr_genes_str"].nunique(dropna=False) == 1 and bool(amr_genes)

        clusters.append({
            "pLIN": plin_code,
            "n_plasmids": len(group),
            "plasmids": plasmids,
            "amr_genes": amr_genes,
            "n_amr_genes": len(amr_genes),
            "amr_profile_consistent": amr_consistent,
            "risk_level": "HIGH" if len(amr_genes) >= 3 else ("MODERATE" if amr_genes else "INFO"),
        })

    # Sort by cluster size first (the primary evidence), then AMR burden
    clusters.sort(key=lambda c: (-c["n_plasmids"], -c["n_amr_genes"]))
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
    cv_data = load_cv_metrics()
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
            contig_count = 0
            for rec in SeqIO.parse(tmp_path, "fasta"):
                seq = str(rec.seq)
                contig_count += 1
                rec_dict = {
                    "plasmid_id": rec.id,
                    "sequence": seq,
                    "length": len(seq),
                    "source_file": uf.name,
                    "_strain_id": uf.name,       # Group contigs from same file
                    "_contig_index": contig_count, # Position within the file
                }
                if auto_detect:
                    result = classify_inc_group(seq, group_names, classifier)
                    rec_dict["inc_type"] = result["predicted_group"]
                    rec_dict["inc_best_match"] = result["best_match"]
                    rec_dict["inc_confidence"] = round(result["confidence"], 4)
                    rec_dict["inc_is_low_confidence"] = result["is_low_confidence"]
                    rec_dict["inc_is_multiple"] = result["is_multiple_inc"]
                    rec_dict["inc_multiple_types"] = result["multiple_inc_types"]
                    rec_dict["inc_top5_candidates"] = result["top5_candidates"]
                    rec_dict["inc_probabilities"] = result["all_probabilities"]
                    rec_dict["inc_blast_used"] = result.get("blast_used", False)
                    rec_dict["inc_blast_hits"] = result.get("blast_replicon_hits", [])
                    # Enrich with CV metrics
                    best = result["best_match"]
                    if best in cv_data["per_class"]:
                        rec_dict["inc_cv_f1"] = cv_data["per_class"][best]["f1"]
                    else:
                        rec_dict["inc_cv_f1"] = None
                    # Check if top-2 are a known confusion pair
                    top2 = result["top5_candidates"][:2]
                    if len(top2) == 2:
                        t2_set = {top2[0][0], top2[1][0]}
                        for cp in cv_data["confusion_pairs"]:
                            if {cp["inc_a"], cp["inc_b"]} == t2_set:
                                rec_dict["inc_confusion_pair"] = True
                                rec_dict["inc_confusion_note"] = (
                                    f"{cp['inc_a']}/{cp['inc_b']} are a known confusion pair "
                                    f"({cp['errors']} CV errors)")
                                break
                        else:
                            rec_dict["inc_confusion_pair"] = False
                            rec_dict["inc_confusion_note"] = ""
                    else:
                        rec_dict["inc_confusion_pair"] = False
                        rec_dict["inc_confusion_note"] = ""
                else:
                    rec_dict["inc_type"] = inc_type

                # Run input quality validation
                rec_dict["quality_warnings"] = validate_sequence_quality(seq, rec.id)

                records.append(rec_dict)
        finally:
            os.unlink(tmp_path)
    return records


@st.cache_data(show_spinner=False)
def compute_kmer_vectors(sequences, k=4):
    """Compute normalised tetranucleotide frequency vectors."""
    bases = "ACGT"
    all_kmers = ["".join(p) for p in iter_product(bases, repeat=k)]
    vectors = np.zeros((len(sequences), len(all_kmers)), dtype=np.float64)
    for idx, seq in enumerate(sequences):
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


# ══════════════════════════════════════════════════════════════════════════════
#  L8: CLUSTER STABILITY ASSESSMENT
# ══════════════════════════════════════════════════════════════════════════════

def assess_cluster_stability(vectors, plin_codes, thresholds=None, n_bootstrap=50):
    """Assess cluster stability via bootstrap resampling.

    Resamples the distance matrix n_bootstrap times, re-clusters, and
    computes co-clustering frequency for each pair of plasmids.

    Returns DataFrame with stability scores per cluster at each level.
    """
    active_thresholds = thresholds if thresholds else PLIN_THRESHOLDS
    n = vectors.shape[0]
    if n < 4:
        return pd.DataFrame()

    rng = np.random.default_rng(42)

    # Track co-clustering at each level
    level_names = list(active_thresholds.keys())
    co_cluster = {level: np.zeros((n, n)) for level in level_names}

    for _ in range(n_bootstrap):
        # Bootstrap resample indices
        boot_idx = rng.choice(n, n, replace=True)
        boot_vectors = vectors[boot_idx]

        boot_dist = pdist(boot_vectors.astype(np.float64), metric="cosine")
        Z_boot = linkage(boot_dist, method="single")

        for level, thresh in active_thresholds.items():
            clusters_boot = fcluster(Z_boot, t=thresh, criterion="distance")
            # Map back to original indices
            for i in range(n):
                for j in range(i + 1, n):
                    orig_i = boot_idx[i]
                    orig_j = boot_idx[j]
                    if clusters_boot[i] == clusters_boot[j]:
                        co_cluster[level][orig_i, orig_j] += 1
                        co_cluster[level][orig_j, orig_i] += 1

    # Compute stability per cluster at each level
    # Parse pLIN codes to get cluster assignments
    plin_parts = [code.split(".") for code in plin_codes]
    results = []
    for li, level in enumerate(level_names):
        cluster_ids = [parts[li] for parts in plin_parts]
        unique_clusters = set(cluster_ids)
        for cl in unique_clusters:
            members = [i for i, c in enumerate(cluster_ids) if c == cl]
            if len(members) < 2:
                results.append({"level": level, "cluster_id": cl,
                                "n_members": 1, "stability_score": 100.0,
                                "stability_status": "Stable"})
                continue

            # Min pairwise co-clustering %
            min_support = 100.0
            for i in range(len(members)):
                for j in range(i + 1, len(members)):
                    support = 100.0 * co_cluster[level][members[i], members[j]] / n_bootstrap
                    min_support = min(min_support, support)

            if min_support >= 80:
                status = "Stable"
            elif min_support >= 50:
                status = "Moderate"
            else:
                status = "Unstable"

            results.append({
                "level": level,
                "cluster_id": cl,
                "n_members": len(members),
                "stability_score": round(min_support, 1),
                "stability_status": status,
            })

    return pd.DataFrame(results) if results else pd.DataFrame()


def compare_linkage_methods(vectors, thresholds=None):
    """Compare clustering results across linkage methods.

    Runs clustering with single, complete, and average linkage and
    computes Adjusted Rand Index between each pair.

    Returns dict with ARI scores and cluster counts per method.
    """
    from sklearn.metrics import adjusted_rand_score
    active_thresholds = thresholds if thresholds else PLIN_THRESHOLDS

    dist_condensed = pdist(vectors.astype(np.float64), metric="cosine")

    methods = ["single", "complete", "average"]
    all_clusters = {}

    for method in methods:
        Z = linkage(dist_condensed, method=method)
        clusters = {}
        for level, thresh in active_thresholds.items():
            clusters[level] = fcluster(Z, t=thresh, criterion="distance")
        all_clusters[method] = clusters

    # Compute ARI between methods at each level
    results = {"level": [], "single_vs_complete": [], "single_vs_average": [],
               "complete_vs_average": [], "n_clusters_single": [],
               "n_clusters_complete": [], "n_clusters_average": []}

    for level in active_thresholds:
        s = all_clusters["single"][level]
        c = all_clusters["complete"][level]
        a = all_clusters["average"][level]

        results["level"].append(level)
        results["single_vs_complete"].append(round(adjusted_rand_score(s, c), 3))
        results["single_vs_average"].append(round(adjusted_rand_score(s, a), 3))
        results["complete_vs_average"].append(round(adjusted_rand_score(c, a), 3))
        results["n_clusters_single"].append(len(set(s)))
        results["n_clusters_complete"].append(len(set(c)))
        results["n_clusters_average"].append(len(set(a)))

    return pd.DataFrame(results)


PLIN_ASSIGNMENTS_PATH = os.path.join(_APP_DIR, "output", "pLIN_assignments.tsv")
REFERENCE_VECTORS_PATH = os.path.join(_APP_DIR, "output", "reference_kmer_vectors.npz")
REFERENCE_ASSIGNMENTS_PATH = os.path.join(_APP_DIR, "output", "pLIN_reference_assignments.tsv")


@st.cache_data(show_spinner=False)
def _load_reference_for_query():
    """Load reference vectors and pLIN codes for nearest-neighbour query mode.

    Uses the full reference database (72,556+ plasmids from PLSDB + training)
    rather than just the smaller training set, so that nearest-neighbour lookups
    match the published case-study results (e.g. Swiss VIM-1 outbreak validation).
    Vectors and pLIN assignments are joined explicitly by plasmid_id, since the
    two files are not stored in the same row order.
    """
    if os.path.exists(REFERENCE_VECTORS_PATH) and os.path.exists(REFERENCE_ASSIGNMENTS_PATH):
        ref_vecs = np.load(REFERENCE_VECTORS_PATH, allow_pickle=True)
        X_ref = ref_vecs["vectors"].astype(np.float64)
        ref_ids = ref_vecs["ids"]

        assign_df = pd.read_csv(REFERENCE_ASSIGNMENTS_PATH, sep="\t")
        assign_df = assign_df.set_index("plasmid_id")

        # Align assignment rows to the vector order by plasmid_id; drop any
        # vectors that have no corresponding pLIN assignment.
        aligned = assign_df.reindex(ref_ids)
        keep_mask = aligned["pLIN"].notna().to_numpy()
        X_ref = X_ref[keep_mask]
        plin_df = aligned.loc[keep_mask].reset_index().rename(columns={"index": "plasmid_id"})
        return X_ref, plin_df

    if not os.path.exists(CLASSIFIER_PATH) or not os.path.exists(PLIN_ASSIGNMENTS_PATH):
        return None, None
    ref_data = np.load(CLASSIFIER_PATH, allow_pickle=True)
    X_ref = ref_data["X"].astype(np.float64)
    plin_df = pd.read_csv(PLIN_ASSIGNMENTS_PATH, sep="\t")
    # Ensure X_ref and plin_df have the same number of rows.
    # After classifier expansion (e.g. Gram-positive groups), X_ref may have
    # more rows than plin_df if assignments haven't been regenerated yet.
    n_ref = min(len(X_ref), len(plin_df))
    X_ref = X_ref[:n_ref]
    plin_df = plin_df.iloc[:n_ref].reset_index(drop=True)
    return X_ref, plin_df


def assign_plin_query_mode(query_vectors, query_records, thresholds=None):
    """Assign pLIN codes to query plasmids by nearest-neighbour lookup.

    Implements the LIN nearest-neighbour assignment rule (Vinatzer et al. 2017):
    for each query, find the nearest reference plasmid and inherit its pLIN code
    at levels where distance <= threshold, creating new branch IDs where it diverges.

    This is mathematically equivalent to single-linkage clustering — not an
    approximation. A query Q belongs to cluster C at threshold t if and only if
    its nearest neighbour N is in C and d(Q,N) <= t.

    Returns:
        plin_codes: list of pLIN code strings
        cluster_assignments: dict {level: array} for compatibility with build_results_df
        query_metadata: list of dicts with nn_plasmid, nn_distance, nn_inc_type
    """
    active_thresholds = thresholds if thresholds else PLIN_THRESHOLDS

    X_ref, plin_df = _load_reference_for_query()
    if X_ref is None:
        raise FileNotFoundError(
            "Reference data not found. Ensure data/inc_classifier.npz and "
            "output/pLIN_assignments.tsv exist for query mode."
        )

    # Max existing cluster IDs at each level (for new branch allocation)
    bin_labels = list(active_thresholds.keys())
    max_ids = {}
    for level in bin_labels:
        max_ids[level] = int(plin_df[f"bin_{level}"].max())

    # Compute cosine distances: each query vs all training references
    dists = cdist(query_vectors, X_ref, metric="cosine")  # (n_query, n_train)

    plin_codes = []
    cluster_assignments = {b: [] for b in bin_labels}
    query_metadata = []

    thresholds_list = list(active_thresholds.values())

    for i in range(len(query_records)):
        nn_idx = np.argmin(dists[i])
        nn_dist = float(dists[i, nn_idx])
        nn_row = plin_df.iloc[nn_idx]

        # Build pLIN code level by level
        code_parts = []
        diverged = False
        for level, thresh in zip(bin_labels, thresholds_list):
            if not diverged and nn_dist <= thresh:
                # Inherit neighbour's code at this level
                bin_val = int(nn_row[f"bin_{level}"])
            else:
                # New branch — assign new unique ID
                diverged = True
                max_ids[level] += 1
                bin_val = max_ids[level]
            code_parts.append(str(bin_val))
            cluster_assignments[level].append(bin_val)

        plin_codes.append(".".join(code_parts))

        query_metadata.append({
            "nn_plasmid": nn_row["plasmid_id"],
            "nn_distance": round(nn_dist, 6),
            "nn_inc_type": nn_row["inc_type"],
            "nn_plin": nn_row["pLIN"],
        })

    # Convert lists to arrays for compatibility
    for level in bin_labels:
        cluster_assignments[level] = np.array(cluster_assignments[level])

    return plin_codes, cluster_assignments, query_metadata


def build_results_df(records, plin_codes, cluster_assignments):
    """Build results DataFrame."""
    # Validate array lengths match
    n_records = len(records)
    n_plin = len(plin_codes)
    if n_records != n_plin:
        raise ValueError(
            f"Mismatch: {n_records} records but {n_plin} pLIN codes. "
            "This may be a caching issue - try refreshing the page."
        )

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
        if "inc_is_multiple" in rec:
            row["inc_is_multiple"] = rec["inc_is_multiple"]
        if "inc_multiple_types" in rec and rec["inc_multiple_types"]:
            # Format multiple Inc types: "IncFII (45%), IncN (30%)"
            multi_str = ", ".join([f"{inc} ({conf*100:.1f}%)" for inc, conf in rec["inc_multiple_types"]])
            row["inc_multiple_types"] = multi_str
        if "inc_best_match" in rec:
            row["inc_best_match"] = rec["inc_best_match"]
        if "inc_top5_candidates" in rec:
            # Format top 5 as string for display: "IncFII (45%), IncN (30%), ..."
            top5_str = ", ".join([f"{inc} ({conf*100:.1f}%)" for inc, conf in rec["inc_top5_candidates"]])
            row["inc_top5_candidates"] = top5_str
        # Track merged contig info
        if "_merged_from" in rec:
            row["merged_contigs"] = ", ".join(rec["_merged_from"])
            row["merged_count"] = rec["_merged_count"]
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
#  MOBSUITE MOBILITY TYPING FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def detect_mobsuite():
    """Auto-detect MOBsuite (mob_typer) binary."""
    binary = None

    # 1. Check PATH
    try:
        result = subprocess.run(["which", "mob_typer"], capture_output=True, text=True)
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
                    candidate = os.path.join(base, env, "bin", "mob_typer")
                    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                        binary = candidate
                        break
            if binary:
                break

    return binary


def run_mobsuite_on_files(uploaded_files, binary, progress_callback=None):
    """Run MOBsuite mob_typer on uploaded FASTA files.

    Returns DataFrame with columns: source_file, predicted_mobility,
    relaxase_type(s), mpf_type, orit_type(s), rep_type(s), etc.
    """
    all_results = []
    total = len(uploaded_files)

    with tempfile.TemporaryDirectory() as tmpdir:
        for idx, uf in enumerate(uploaded_files):
            fasta_path = os.path.join(tmpdir, uf.name)
            with open(fasta_path, "wb") as f:
                f.write(uf.getvalue())

            out_path = os.path.join(tmpdir, f"{uf.name}.mobtyper.txt")
            source_name = uf.name.replace(".fasta", "").replace(".fa", "").replace(".fna", "")

            cmd = [binary, "--infile", fasta_path, "--out_file", out_path]

            try:
                subprocess.run(cmd, capture_output=True, timeout=300)
                if os.path.isfile(out_path):
                    df = pd.read_csv(out_path, sep="\t")
                    df.insert(0, "source_file", source_name)
                    all_results.append(df)
            except Exception:
                pass

            if progress_callback:
                progress_callback((idx + 1) / total)

    if all_results:
        return pd.concat(all_results, ignore_index=True)
    return pd.DataFrame()


# ══════════════════════════════════════════════════════════════════════════════
#  MASH / MINHASH ANI ESTIMATION
# ══════════════════════════════════════════════════════════════════════════════

def detect_mash():
    """Auto-detect Mash binary."""
    binary = None
    try:
        result = subprocess.run(["which", "mash"], capture_output=True, text=True)
        if result.returncode == 0 and result.stdout.strip():
            binary = result.stdout.strip()
    except Exception:
        pass
    if not binary:
        home = os.path.expanduser("~")
        for base in [
            os.path.join(home, "miniconda3", "envs"),
            os.path.join(home, "miniforge3", "envs"),
            os.path.join(home, "anaconda3", "envs"),
            os.path.join(home, "mambaforge", "envs"),
        ]:
            if os.path.isdir(base):
                for env in sorted(os.listdir(base)):
                    candidate = os.path.join(base, env, "bin", "mash")
                    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                        binary = candidate
                        break
            if binary:
                break
    return binary


def run_mash_distances(uploaded_files, mash_binary, progress_callback=None):
    """Run Mash pairwise distance estimation on uploaded FASTA files.

    Returns DataFrame with columns: query, reference, mash_distance, p_value, matching_hashes, ani_estimate
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write all FASTAs to temp dir
        fasta_paths = []
        for uf in uploaded_files:
            fasta_path = os.path.join(tmpdir, uf.name)
            with open(fasta_path, "wb") as f:
                f.write(uf.getvalue())
            fasta_paths.append(fasta_path)

        if len(fasta_paths) < 2:
            return pd.DataFrame()

        # Create combined sketch
        sketch_path = os.path.join(tmpdir, "all_plasmids.msh")
        cmd_sketch = [mash_binary, "sketch", "-o", sketch_path, "-k", "21", "-s", "10000"] + fasta_paths
        try:
            subprocess.run(cmd_sketch, capture_output=True, timeout=600)
        except Exception:
            return pd.DataFrame()

        # Pairwise distances
        dist_path = os.path.join(tmpdir, "mash_dist.tsv")
        cmd_dist = [mash_binary, "dist", sketch_path, sketch_path]
        try:
            result = subprocess.run(cmd_dist, capture_output=True, text=True, timeout=600)
            if result.returncode != 0:
                return pd.DataFrame()
        except Exception:
            return pd.DataFrame()

        # Parse Mash output: ref\tquery\tdist\tp-value\tmatching-hashes
        rows = []
        for line in result.stdout.strip().split("\n"):
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 5:
                ref = os.path.basename(parts[0]).replace(".fasta", "").replace(".fa", "").replace(".fna", "")
                query = os.path.basename(parts[1]).replace(".fasta", "").replace(".fa", "").replace(".fna", "")
                if ref == query:
                    continue
                mash_dist = float(parts[2])
                p_value = float(parts[3])
                matching = parts[4]
                ani_est = round((1 - mash_dist) * 100, 2)
                rows.append({
                    "query": query, "reference": ref,
                    "mash_distance": round(mash_dist, 6),
                    "p_value": p_value,
                    "matching_hashes": matching,
                    "ani_estimate": ani_est,
                })

        if progress_callback:
            progress_callback(1.0)

    return pd.DataFrame(rows) if rows else pd.DataFrame()


def check_ani_concordance(plin_df, mash_df):
    """Cross-validate pLIN cosine NN distances against Mash ANI estimates.

    For each plasmid pair in the Mash results, check if cosine distance
    and Mash ANI agree. Discordant cases suggest the 4-mer composition
    proxy may be unreliable for that pair.

    Returns DataFrame with concordance status per plasmid.
    """
    if mash_df is None or len(mash_df) == 0 or plin_df is None:
        return pd.DataFrame()

    # Build a lookup of NN distance per plasmid from plin_df
    nn_lookup = {}
    if "nn_distance" in plin_df.columns and "nn_plasmid" in plin_df.columns:
        for _, row in plin_df.iterrows():
            nn_lookup[row["plasmid_id"]] = {
                "nn_distance": row.get("nn_distance", None),
                "nn_plasmid": row.get("nn_plasmid", ""),
            }

    # Build Mash ANI lookup for pairs
    mash_pairs = {}
    for _, row in mash_df.iterrows():
        q, r = str(row["query"]), str(row["reference"])
        ani = row.get("ani_estimate", 0)
        # Store both directions
        mash_pairs[(q, r)] = ani
        mash_pairs[(r, q)] = ani

    results = []
    for pid, nn_info in nn_lookup.items():
        nn_dist = nn_info["nn_distance"]
        nn_id = nn_info["nn_plasmid"]
        if nn_dist is None:
            continue

        # Try to find this pair in Mash results
        mash_ani = mash_pairs.get((pid, nn_id))
        if mash_ani is None:
            # Try partial ID match (Mash may use filename-based IDs)
            for (q, r), ani in mash_pairs.items():
                if pid in q and nn_id in r:
                    mash_ani = ani
                    break

        if mash_ani is not None:
            # Concordance logic:
            # Low distance + low ANI = discordant (composition similar but sequence divergent)
            # High distance + high ANI = discordant (rare, but possible with rearrangements)
            if nn_dist < 0.01 and mash_ani < 95:
                status = "discordant"
                note = f"Low cosine distance ({nn_dist:.4f}) but low ANI ({mash_ani:.1f}%)"
            elif nn_dist > 0.02 and mash_ani > 98:
                status = "discordant"
                note = f"High cosine distance ({nn_dist:.4f}) but high ANI ({mash_ani:.1f}%)"
            else:
                status = "concordant"
                note = ""
        else:
            status = "no_mash_data"
            note = "No Mash pair data for this NN"

        results.append({
            "plasmid_id": pid,
            "nn_distance": round(nn_dist, 6) if nn_dist else None,
            "mash_ani": round(mash_ani, 1) if mash_ani else None,
            "concordance": status,
            "note": note,
        })

    return pd.DataFrame(results) if results else pd.DataFrame()


# ══════════════════════════════════════════════════════════════════════════════
#  FASTANI INTEGRATION
# ══════════════════════════════════════════════════════════════════════════════

def detect_fastani():
    """Auto-detect FastANI binary."""
    binary = None
    try:
        result = subprocess.run(["which", "fastANI"], capture_output=True, text=True)
        if result.returncode == 0 and result.stdout.strip():
            binary = result.stdout.strip()
    except Exception:
        pass
    if not binary:
        home = os.path.expanduser("~")
        for base in [
            os.path.join(home, "miniconda3", "envs"),
            os.path.join(home, "miniforge3", "envs"),
            os.path.join(home, "anaconda3", "envs"),
            os.path.join(home, "mambaforge", "envs"),
        ]:
            if os.path.isdir(base):
                for env in sorted(os.listdir(base)):
                    candidate = os.path.join(base, env, "bin", "fastANI")
                    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                        binary = candidate
                        break
            if binary:
                break
    return binary


def run_fastani(uploaded_files, fastani_binary, progress_callback=None):
    """Run FastANI all-vs-all on uploaded FASTA files.

    Returns DataFrame with columns: query, reference, ani, orthologous_matches, total_fragments
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_paths = []
        for uf in uploaded_files:
            fasta_path = os.path.join(tmpdir, uf.name)
            with open(fasta_path, "wb") as f:
                f.write(uf.getvalue())
            fasta_paths.append(fasta_path)

        if len(fasta_paths) < 2:
            return pd.DataFrame()

        # Write query and reference lists
        list_path = os.path.join(tmpdir, "file_list.txt")
        with open(list_path, "w") as f:
            for p in fasta_paths:
                f.write(p + "\n")

        out_path = os.path.join(tmpdir, "fastani_out.tsv")
        cmd = [
            fastani_binary,
            "--ql", list_path,
            "--rl", list_path,
            "-o", out_path,
            "--fragLen", "1000",  # Smaller fragment for plasmids
            "-t", "4",
        ]

        try:
            subprocess.run(cmd, capture_output=True, timeout=1800)
        except Exception:
            return pd.DataFrame()

        if not os.path.isfile(out_path) or os.path.getsize(out_path) == 0:
            return pd.DataFrame()

        # Parse FastANI output: query\treference\tANI\torthologous_matches\ttotal_fragments
        rows = []
        with open(out_path, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    query = os.path.basename(parts[0]).replace(".fasta", "").replace(".fa", "").replace(".fna", "")
                    ref = os.path.basename(parts[1]).replace(".fasta", "").replace(".fa", "").replace(".fna", "")
                    if query == ref:
                        continue
                    rows.append({
                        "query": query, "reference": ref,
                        "ani": round(float(parts[2]), 2),
                        "orthologous_matches": int(parts[3]),
                        "total_fragments": int(parts[4]),
                    })

        if progress_callback:
            progress_callback(1.0)

    return pd.DataFrame(rows) if rows else pd.DataFrame()


# ══════════════════════════════════════════════════════════════════════════════
#  SNP SUB-TYPING WITHIN L6 CLUSTERS
# ══════════════════════════════════════════════════════════════════════════════

def detect_minimap2():
    """Auto-detect minimap2 binary."""
    binary = None
    try:
        result = subprocess.run(["which", "minimap2"], capture_output=True, text=True)
        if result.returncode == 0 and result.stdout.strip():
            binary = result.stdout.strip()
    except Exception:
        pass
    if not binary:
        home = os.path.expanduser("~")
        for base in [
            os.path.join(home, "miniconda3", "envs"),
            os.path.join(home, "miniforge3", "envs"),
            os.path.join(home, "anaconda3", "envs"),
            os.path.join(home, "mambaforge", "envs"),
        ]:
            if os.path.isdir(base):
                for env in sorted(os.listdir(base)):
                    candidate = os.path.join(base, env, "bin", "minimap2")
                    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                        binary = candidate
                        break
            if binary:
                break
    return binary


def run_snp_subtyping(records, cluster_assignments, minimap2_binary, progress_callback=None):
    """Run SNP-level sub-typing within each L6 cluster using minimap2 + paftools variant calling.

    For each L6 cluster with >=2 plasmids, aligns all members against the longest
    member as reference and counts mismatches (approximate SNP differences).

    Returns DataFrame with columns: l6_cluster, plasmid_id, reference_id, snp_count, alignment_identity
    """
    l6_clusters = {}
    for i, rec in enumerate(records):
        cl = int(cluster_assignments["F"][i])
        l6_clusters.setdefault(cl, []).append(i)

    # Only process clusters with >=2 members
    multi_clusters = {k: v for k, v in l6_clusters.items() if len(v) >= 2}
    if not multi_clusters:
        return pd.DataFrame()

    all_results = []
    total = len(multi_clusters)

    with tempfile.TemporaryDirectory() as tmpdir:
        for ci, (cluster_id, member_indices) in enumerate(multi_clusters.items()):
            # Pick longest plasmid as reference
            ref_idx = max(member_indices, key=lambda i: records[i]["length"])
            ref_rec = records[ref_idx]

            ref_path = os.path.join(tmpdir, f"ref_c{cluster_id}.fasta")
            with open(ref_path, "w") as f:
                f.write(f">{ref_rec['plasmid_id']}\n{ref_rec['sequence']}\n")

            for m_idx in member_indices:
                if m_idx == ref_idx:
                    all_results.append({
                        "l6_cluster": cluster_id,
                        "plasmid_id": records[m_idx]["plasmid_id"],
                        "reference_id": ref_rec["plasmid_id"],
                        "snp_count": 0,
                        "alignment_identity": 100.0,
                    })
                    continue

                query_rec = records[m_idx]
                query_path = os.path.join(tmpdir, f"q_c{cluster_id}_{m_idx}.fasta")
                with open(query_path, "w") as f:
                    f.write(f">{query_rec['plasmid_id']}\n{query_rec['sequence']}\n")

                # Run minimap2 -cx asm5 (for closely related sequences)
                try:
                    result = subprocess.run(
                        [minimap2_binary, "-cx", "asm5", "--cs", ref_path, query_path],
                        capture_output=True, text=True, timeout=120,
                    )
                    if result.returncode == 0 and result.stdout.strip():
                        # Parse PAF output to extract alignment identity and mismatches
                        best_identity = 0.0
                        total_mismatches = 0
                        total_aligned = 0
                        for line in result.stdout.strip().split("\n"):
                            fields = line.split("\t")
                            if len(fields) >= 12:
                                matches = int(fields[9])
                                block_len = int(fields[10])
                                mismatches = block_len - matches
                                identity = matches / block_len * 100 if block_len > 0 else 0
                                if block_len > total_aligned:
                                    best_identity = identity
                                    total_mismatches = mismatches
                                    total_aligned = block_len

                        all_results.append({
                            "l6_cluster": cluster_id,
                            "plasmid_id": query_rec["plasmid_id"],
                            "reference_id": ref_rec["plasmid_id"],
                            "snp_count": total_mismatches,
                            "alignment_identity": round(best_identity, 3),
                        })
                    else:
                        all_results.append({
                            "l6_cluster": cluster_id,
                            "plasmid_id": query_rec["plasmid_id"],
                            "reference_id": ref_rec["plasmid_id"],
                            "snp_count": -1,
                            "alignment_identity": 0.0,
                        })
                except Exception:
                    all_results.append({
                        "l6_cluster": cluster_id,
                        "plasmid_id": query_rec["plasmid_id"],
                        "reference_id": ref_rec["plasmid_id"],
                        "snp_count": -1,
                        "alignment_identity": 0.0,
                    })

            if progress_callback:
                progress_callback((ci + 1) / total)

    return pd.DataFrame(all_results) if all_results else pd.DataFrame()


# ══════════════════════════════════════════════════════════════════════════════
#  L5: RECOMBINATION DETECTION
# ══════════════════════════════════════════════════════════════════════════════

def detect_recombination_signals(records, cluster_assignments, minimap2_binary,
                                  progress_callback=None):
    """Detect potential recombination events by analyzing alignment coverage patterns.

    For each L6 cluster with >=2 members, runs minimap2 and analyzes:
    - Alignment coverage (what fraction of query aligns to reference)
    - Number of alignment blocks (fragmented = possible rearrangement)
    - Largest unaligned region (putative recombination breakpoint)

    Returns DataFrame with recombination flags per plasmid.
    """
    l6_clusters = {}
    for i, rec in enumerate(records):
        cl = int(cluster_assignments["F"][i])
        l6_clusters.setdefault(cl, []).append(i)

    multi_clusters = {k: v for k, v in l6_clusters.items() if len(v) >= 2}
    if not multi_clusters:
        return pd.DataFrame()

    all_results = []
    total = len(multi_clusters)

    with tempfile.TemporaryDirectory() as tmpdir:
        for ci, (cluster_id, member_indices) in enumerate(multi_clusters.items()):
            ref_idx = max(member_indices, key=lambda i: records[i]["length"])
            ref_rec = records[ref_idx]
            ref_len = ref_rec["length"]

            ref_path = os.path.join(tmpdir, f"ref_c{cluster_id}.fasta")
            with open(ref_path, "w") as f:
                f.write(f">{ref_rec['plasmid_id']}\n{ref_rec['sequence']}\n")

            for m_idx in member_indices:
                if m_idx == ref_idx:
                    all_results.append({
                        "plasmid_id": ref_rec["plasmid_id"],
                        "l6_cluster": cluster_id,
                        "alignment_coverage_pct": 100.0,
                        "n_alignment_blocks": 1,
                        "largest_unaligned_bp": 0,
                        "recombination_flag": "None",
                    })
                    continue

                query_rec = records[m_idx]
                query_len = query_rec["length"]
                query_path = os.path.join(tmpdir, f"q_c{cluster_id}_{m_idx}.fasta")
                with open(query_path, "w") as f:
                    f.write(f">{query_rec['plasmid_id']}\n{query_rec['sequence']}\n")

                try:
                    result = subprocess.run(
                        [minimap2_binary, "-cx", "asm5", ref_path, query_path],
                        capture_output=True, text=True, timeout=120,
                    )
                    if result.returncode == 0 and result.stdout.strip():
                        # Parse PAF: collect alignment blocks
                        blocks = []
                        for line in result.stdout.strip().split("\n"):
                            fields = line.split("\t")
                            if len(fields) >= 12:
                                q_start = int(fields[2])
                                q_end = int(fields[3])
                                matches = int(fields[9])
                                block_len = int(fields[10])
                                blocks.append((q_start, q_end, matches, block_len))

                        if blocks:
                            # Merge overlapping blocks
                            blocks.sort()
                            merged = [blocks[0]]
                            for b in blocks[1:]:
                                if b[0] <= merged[-1][1]:
                                    merged[-1] = (merged[-1][0], max(merged[-1][1], b[1]),
                                                  merged[-1][2] + b[2], merged[-1][3] + b[3])
                                else:
                                    merged.append(b)

                            total_aligned = sum(b[1] - b[0] for b in merged)
                            coverage = 100.0 * total_aligned / query_len if query_len > 0 else 0

                            # Find gaps between blocks
                            gaps = []
                            for k in range(1, len(merged)):
                                gap = merged[k][0] - merged[k - 1][1]
                                if gap > 100:
                                    gaps.append(gap)
                            largest_gap = max(gaps) if gaps else 0

                            # Recombination flag
                            if coverage < 50:
                                flag = "High"
                            elif coverage < 70 or largest_gap > 5000:
                                flag = "Medium"
                            elif len(merged) > 3 or largest_gap > 2000:
                                flag = "Low"
                            else:
                                flag = "None"

                            all_results.append({
                                "plasmid_id": query_rec["plasmid_id"],
                                "l6_cluster": cluster_id,
                                "alignment_coverage_pct": round(coverage, 1),
                                "n_alignment_blocks": len(merged),
                                "largest_unaligned_bp": largest_gap,
                                "recombination_flag": flag,
                            })
                        else:
                            all_results.append({
                                "plasmid_id": query_rec["plasmid_id"],
                                "l6_cluster": cluster_id,
                                "alignment_coverage_pct": 0.0,
                                "n_alignment_blocks": 0,
                                "largest_unaligned_bp": query_len,
                                "recombination_flag": "High",
                            })
                    else:
                        all_results.append({
                            "plasmid_id": query_rec["plasmid_id"],
                            "l6_cluster": cluster_id,
                            "alignment_coverage_pct": 0.0,
                            "n_alignment_blocks": 0,
                            "largest_unaligned_bp": query_len,
                            "recombination_flag": "High",
                        })
                except Exception:
                    all_results.append({
                        "plasmid_id": query_rec["plasmid_id"],
                        "l6_cluster": cluster_id,
                        "alignment_coverage_pct": 0.0,
                        "n_alignment_blocks": 0,
                        "largest_unaligned_bp": 0,
                        "recombination_flag": "Unknown",
                    })

            if progress_callback:
                progress_callback((ci + 1) / total)

    return pd.DataFrame(all_results) if all_results else pd.DataFrame()


# ══════════════════════════════════════════════════════════════════════════════
#  TEMPORAL OUTBREAK CLUSTERING
# ══════════════════════════════════════════════════════════════════════════════

def detect_temporal_outbreak_clusters(plin_df, integrated_df, metadata_df, time_window_days=30):
    """Flag potential outbreak clusters with temporal evidence.

    Extends basic outbreak detection by requiring plasmids to share:
    1. Same L6 pLIN cluster (bin_F)
    2. Identical AMR resistance profile
    3. Collection dates within a specified time window

    Returns list of dicts with cluster info including temporal evidence.
    """
    if integrated_df is None or len(integrated_df) == 0 or metadata_df is None:
        return []

    df = integrated_df.copy()

    # Merge metadata date columns
    date_col = None
    for col in metadata_df.columns:
        if any(kw in col.lower() for kw in ["date", "collection_date", "sample_date"]):
            date_col = col
            break

    if date_col is None:
        return []

    meta_dates = metadata_df[["plasmid_id", date_col]].copy()
    meta_dates = meta_dates.rename(columns={date_col: "_collection_date"})
    meta_dates["_collection_date"] = pd.to_datetime(meta_dates["_collection_date"], errors="coerce")
    df = df.merge(meta_dates, on="plasmid_id", how="left")

    if df["_collection_date"].isna().all():
        return []

    # Gather AMR genes per plasmid, purely as supplementary evidence (not required).
    if "AMR_genes" in df.columns:
        df["_amr_genes_str"] = df["AMR_genes"].fillna("")
    else:
        df["_amr_genes_str"] = ""

    # Merge location if available
    location_col = None
    for col in metadata_df.columns:
        if any(kw in col.lower() for kw in ["location", "ward", "hospital", "site", "unit"]):
            location_col = col
            break

    clusters = []
    for plin_code, group in df.groupby("pLIN"):
        if len(group) < 2 or not isinstance(plin_code, str):
            continue
        if plin_code == "UNCLASSIFIED" or "NEW" in plin_code:
            continue  # Divergent/novel codes are not a shared-plasmid match

        # Filter to those with valid dates
        dated = group.dropna(subset=["_collection_date"])
        if len(dated) < 2:
            continue

        # Check if collection dates fall within time window
        date_range = (dated["_collection_date"].max() - dated["_collection_date"].min()).days
        if date_range <= time_window_days:
            plasmids = dated["plasmid_id"].tolist()
            amr_genes = sorted({
                g.strip() for x in dated["_amr_genes_str"] for g in x.split(";") if g.strip()
            })

            # Gather location info
            locations = []
            if location_col and location_col in metadata_df.columns:
                loc_data = dated.merge(
                    metadata_df[["plasmid_id", location_col]], on="plasmid_id", how="left"
                )
                locations = loc_data[location_col].dropna().unique().tolist()

            risk = "CRITICAL" if len(amr_genes) >= 3 and date_range <= 7 else \
                   "HIGH" if len(amr_genes) >= 3 or date_range <= 7 else "MODERATE"

            clusters.append({
                "pLIN": plin_code,
                "n_plasmids": len(dated),
                "plasmids": plasmids,
                "amr_genes": amr_genes,
                "n_amr_genes": len(amr_genes),
                "date_range_days": date_range,
                "earliest_date": str(dated["_collection_date"].min().date()),
                "latest_date": str(dated["_collection_date"].max().date()),
                "locations": locations,
                "risk_level": risk,
                "temporal_evidence": True,
            })

    clusters.sort(key=lambda c: (-{"CRITICAL": 3, "HIGH": 2, "MODERATE": 1}.get(c["risk_level"], 0),
                                  -c["n_amr_genes"], -c["n_plasmids"]))
    return clusters


# ══════════════════════════════════════════════════════════════════════════════
#  L7: EVOLUTIONARY RATE ESTIMATION (MOLECULAR CLOCK)
# ══════════════════════════════════════════════════════════════════════════════

def estimate_evolutionary_rate(snp_df, metadata_df, records):
    """Estimate SNP accumulation rate (molecular clock) for L6 clusters.

    Requires both SNP data (from run_snp_subtyping) and collection dates
    from metadata. For each L6 cluster with >=3 dated members, fits a
    linear regression of SNP count vs time difference.

    Returns DataFrame with per-cluster rate estimates.
    """
    if snp_df is None or len(snp_df) == 0 or metadata_df is None:
        return pd.DataFrame()

    # Find date column
    date_col = None
    for col in metadata_df.columns:
        if any(kw in col.lower() for kw in ["date", "collection_date", "sample_date"]):
            date_col = col
            break
    if date_col is None:
        return pd.DataFrame()

    # Build date lookup
    date_lookup = {}
    for _, row in metadata_df.iterrows():
        pid = str(row.get("plasmid_id", ""))
        try:
            dt = pd.to_datetime(row[date_col])
            if pd.notna(dt):
                date_lookup[pid] = dt
        except Exception:
            pass

    if len(date_lookup) < 3:
        return pd.DataFrame()

    # Build plasmid length lookup
    len_lookup = {r["plasmid_id"]: r["length"] for r in records}

    # Group SNP data by L6 cluster
    results = []
    for cluster_id, group in snp_df.groupby("l6_cluster"):
        # Get dated members
        dated_members = []
        for _, row in group.iterrows():
            pid = row["plasmid_id"]
            if pid in date_lookup and row["snp_count"] >= 0:
                dated_members.append({
                    "plasmid_id": pid,
                    "snp_count": int(row["snp_count"]),
                    "date": date_lookup[pid],
                })

        if len(dated_members) < 3:
            continue

        # Compute pairwise time differences and SNP differences
        from scipy.stats import linregress
        time_diffs = []
        snp_diffs = []
        for i in range(len(dated_members)):
            for j in range(i + 1, len(dated_members)):
                days = abs((dated_members[i]["date"] - dated_members[j]["date"]).days)
                snps = abs(dated_members[i]["snp_count"] - dated_members[j]["snp_count"])
                if days > 0:
                    time_diffs.append(days)
                    snp_diffs.append(snps)

        if len(time_diffs) < 3:
            continue

        # Linear regression
        slope, intercept, r_value, p_value, std_err = linregress(time_diffs, snp_diffs)
        snps_per_year = slope * 365.25

        # Normalize by plasmid length
        ref_pid = dated_members[0]["plasmid_id"]
        plasmid_len = len_lookup.get(ref_pid, 100000)
        subs_per_site_per_year = snps_per_year / plasmid_len if plasmid_len > 0 else 0

        # Interpretation
        if subs_per_site_per_year > 1e-4:
            interpretation = "Unusually high (possible recombination)"
        elif subs_per_site_per_year > 1e-5:
            interpretation = "Within expected range for plasmids"
        elif subs_per_site_per_year > 1e-7:
            interpretation = "Low (conserved backbone)"
        else:
            interpretation = "Very low (possible same-source)"

        results.append({
            "l6_cluster": int(cluster_id),
            "n_dated_members": len(dated_members),
            "n_pairwise_comparisons": len(time_diffs),
            "snps_per_year": round(snps_per_year, 2),
            "subs_per_site_per_year": f"{subs_per_site_per_year:.2e}",
            "r_squared": round(r_value ** 2, 3),
            "p_value": round(p_value, 4),
            "interpretation": interpretation,
        })

    return pd.DataFrame(results) if results else pd.DataFrame()


# ══════════════════════════════════════════════════════════════════════════════
#  CRISPR HOST INFERENCE FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def detect_minced():
    """Auto-detect MinCED (Mining CRISPRs in Environmental Datasets) binary."""
    binary = None
    try:
        result = subprocess.run(["which", "minced"], capture_output=True, text=True)
        if result.returncode == 0 and result.stdout.strip():
            binary = result.stdout.strip()
    except Exception:
        pass

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
                    candidate = os.path.join(base, env, "bin", "minced")
                    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                        binary = candidate
                        break
            if binary:
                break

    if not binary:
        home = os.path.expanduser("~")
        for prefix in [
            os.path.join(home, "miniforge3", "bin", "minced"),
            os.path.join(home, "miniconda3", "bin", "minced"),
            os.path.join(home, "anaconda3", "bin", "minced"),
            os.path.join(home, "mambaforge", "bin", "minced"),
        ]:
            if os.path.isfile(prefix) and os.access(prefix, os.X_OK):
                binary = prefix
                break

    return binary


def detect_blastn():
    """Auto-detect BLAST+ blastn and makeblastdb binaries.

    Returns (blastn_path, makeblastdb_path) — either may be None.
    """
    blastn = None
    makeblastdb = None

    for tool_name in ["blastn", "makeblastdb"]:
        try:
            result = subprocess.run(["which", tool_name], capture_output=True, text=True)
            if result.returncode == 0 and result.stdout.strip():
                if tool_name == "blastn":
                    blastn = result.stdout.strip()
                else:
                    makeblastdb = result.stdout.strip()
        except Exception:
            pass

    home = os.path.expanduser("~")
    search_bases = [
        os.path.join(home, "miniforge3"),
        os.path.join(home, "miniconda3"),
        os.path.join(home, "anaconda3"),
        os.path.join(home, "mambaforge"),
    ]

    for tool_name in ["blastn", "makeblastdb"]:
        current = blastn if tool_name == "blastn" else makeblastdb
        if current:
            continue
        for base in search_bases:
            candidate = os.path.join(base, "bin", tool_name)
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                if tool_name == "blastn":
                    blastn = candidate
                else:
                    makeblastdb = candidate
                break
        current = blastn if tool_name == "blastn" else makeblastdb
        if current:
            continue
        for base in search_bases:
            envs_dir = os.path.join(base, "envs")
            if os.path.isdir(envs_dir):
                for env in sorted(os.listdir(envs_dir)):
                    candidate = os.path.join(envs_dir, env, "bin", tool_name)
                    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                        if tool_name == "blastn":
                            blastn = candidate
                        else:
                            makeblastdb = candidate
                        break
                if (tool_name == "blastn" and blastn) or (tool_name == "makeblastdb" and makeblastdb):
                    break

    return blastn, makeblastdb


# ── MLST chromosomal typing ──────────────────────────────────────────────────


def detect_mlst():
    """Auto-detect mlst (Torsten Seemann's MLST tool) binary."""
    binary = None
    try:
        result = subprocess.run(["which", "mlst"], capture_output=True, text=True)
        if result.returncode == 0 and result.stdout.strip():
            binary = result.stdout.strip()
    except Exception:
        pass

    if not binary:
        home = os.path.expanduser("~")
        for base in [
            os.path.join(home, "miniconda3", "envs"),
            os.path.join(home, "miniforge3", "envs"),
            os.path.join(home, "anaconda3", "envs"),
            os.path.join(home, "mambaforge", "envs"),
        ]:
            if os.path.isdir(base):
                for env in sorted(os.listdir(base)):
                    candidate = os.path.join(base, env, "bin", "mlst")
                    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                        binary = candidate
                        break
            if binary:
                break

    if not binary:
        home = os.path.expanduser("~")
        for prefix in [
            os.path.join(home, "miniforge3", "bin", "mlst"),
            os.path.join(home, "miniconda3", "bin", "mlst"),
        ]:
            if os.path.isfile(prefix) and os.access(prefix, os.X_OK):
                binary = prefix
                break

    return binary


def _mlst_perl_env(binary):
    """Build environment dict with correct PERL5LIB for mlst."""
    env = os.environ.copy()
    # mlst needs Perl libs from its conda env
    bin_dir = os.path.dirname(binary)
    env_root = os.path.dirname(bin_dir)
    perl_lib = os.path.join(env_root, "lib", "perl5", "site_perl")
    perl_lib2 = os.path.join(env_root, "lib", "perl5")
    existing = env.get("PERL5LIB", "")
    env["PERL5LIB"] = f"{perl_lib}:{perl_lib2}:{existing}" if existing else f"{perl_lib}:{perl_lib2}"
    return env


def run_mlst_on_genomes(genome_files, mlst_binary, progress_callback=None):
    """Run MLST typing on host bacterial genome FASTA files.

    Returns DataFrame with columns: genome_name, scheme, ST, alleles.
    """
    all_results = []
    total = len(genome_files)
    env = _mlst_perl_env(mlst_binary)

    with tempfile.TemporaryDirectory() as tmpdir:
        for idx, uf in enumerate(genome_files):
            genome_name = uf.name
            for ext in [".fasta", ".fa", ".fna"]:
                genome_name = genome_name.replace(ext, "")
            fasta_path = os.path.join(tmpdir, uf.name)
            with open(fasta_path, "wb") as f:
                f.write(uf.getvalue())

            try:
                result = subprocess.run(
                    [mlst_binary, fasta_path],
                    capture_output=True, text=True, timeout=300, env=env,
                )
                if result.returncode == 0 and result.stdout.strip():
                    parts = result.stdout.strip().split("\t")
                    scheme = parts[1] if len(parts) > 1 else "-"
                    st_type = parts[2] if len(parts) > 2 else "-"
                    alleles = ";".join(parts[3:]) if len(parts) > 3 else ""
                    all_results.append({
                        "genome_name": genome_name,
                        "scheme": scheme,
                        "ST": st_type,
                        "alleles": alleles,
                    })
                else:
                    all_results.append({
                        "genome_name": genome_name,
                        "scheme": "-",
                        "ST": "-",
                        "alleles": "",
                    })
            except Exception:
                all_results.append({
                    "genome_name": genome_name,
                    "scheme": "-",
                    "ST": "-",
                    "alleles": "",
                })

            if progress_callback:
                progress_callback((idx + 1) / total)

    return pd.DataFrame(all_results) if all_results else pd.DataFrame()


def build_genome_plasmid_mapping(genome_files, plasmid_ids, metadata_df=None):
    """Link genome files to plasmid IDs via filename prefix matching or metadata.

    Strategy: (1) shared filename tokens, (2) metadata CSV with genome column.
    Returns DataFrame with genome_name, plasmid_id columns.
    """
    mappings = []

    def _clean(name):
        for ext in [".fasta", ".fa", ".fna"]:
            name = name.replace(ext, "")
        return name

    genome_names = [_clean(gf.name) for gf in genome_files]
    plasmid_clean = {_clean(pid): pid for pid in plasmid_ids}

    # Strategy 1: shared filename tokens
    for gn in genome_names:
        gn_parts = set(gn.replace("-", "_").split("_")) - {
            "genome", "chromosome", "chr", "contig", "assembly",
        }
        for pc, pid_orig in plasmid_clean.items():
            pc_parts = set(pc.replace("-", "_").split("_")) - {
                "plasmid", "plas", "contig",
            }
            shared = gn_parts & pc_parts
            if shared and len(shared) >= 1:
                mappings.append({"genome_name": gn, "plasmid_id": pid_orig})

    # Strategy 2: metadata-based
    if not mappings and metadata_df is not None:
        genome_col = None
        for col in metadata_df.columns:
            if any(kw in col.lower() for kw in ["genome", "chromosome", "host_genome"]):
                genome_col = col
                break
        if genome_col and "plasmid_id" in metadata_df.columns:
            for _, row in metadata_df.iterrows():
                gname = _clean(str(row[genome_col]))
                pid = str(row["plasmid_id"])
                if gname in genome_names:
                    mappings.append({"genome_name": gname, "plasmid_id": pid})

    return pd.DataFrame(mappings) if mappings else pd.DataFrame(columns=["genome_name", "plasmid_id"])


def classify_transmission_mode(mlst_df, plin_df, genome_plasmid_mapping):
    """Classify pairwise transmission as clonal spread, HGT, or independent.

    Returns (pairs_df, summary_dict).
    """
    merged = genome_plasmid_mapping.merge(mlst_df, on="genome_name", how="left")
    plin_cols = ["plasmid_id"]
    if "pLIN" in plin_df.columns:
        plin_cols.append("pLIN")
    if "bin_F" in plin_df.columns:
        plin_cols.append("bin_F")
    if "inc_type" in plin_df.columns:
        plin_cols.append("inc_type")
    merged = merged.merge(plin_df[plin_cols], on="plasmid_id", how="left")

    pairs = []
    samples = merged.to_dict("records")
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            a, b = samples[i], samples[j]
            st_a = str(a.get("ST", "-"))
            st_b = str(b.get("ST", "-"))
            same_st = st_a == st_b and st_a not in ["-", "", "0", "None"]
            bin_a = a.get("bin_F", "")
            bin_b = b.get("bin_F", "")
            same_plin = bin_a == bin_b and bin_a not in ["", None]

            if same_st and same_plin:
                mode = "Clonal spread"
            elif not same_st and same_plin:
                mode = "Horizontal plasmid transfer"
            elif same_st and not same_plin:
                mode = "Same strain, different plasmids"
            else:
                mode = "Independent"

            pairs.append({
                "genome_A": a.get("genome_name", ""),
                "genome_B": b.get("genome_name", ""),
                "ST_A": st_a,
                "ST_B": st_b,
                "plasmid_A": a.get("plasmid_id", ""),
                "plasmid_B": b.get("plasmid_id", ""),
                "pLIN_A": a.get("pLIN", ""),
                "pLIN_B": b.get("pLIN", ""),
                "same_ST": same_st,
                "same_pLIN_L6": same_plin,
                "transmission_mode": mode,
            })

    pairs_df = pd.DataFrame(pairs) if pairs else pd.DataFrame()
    summary = {}
    if len(pairs_df) > 0:
        mc = pairs_df["transmission_mode"].value_counts()
        summary = {
            "total_pairs": len(pairs_df),
            "clonal_spread": int(mc.get("Clonal spread", 0)),
            "horizontal_transfer": int(mc.get("Horizontal plasmid transfer", 0)),
            "same_strain_diff_plasmid": int(mc.get("Same strain, different plasmids", 0)),
            "independent": int(mc.get("Independent", 0)),
        }
    return pairs_df, summary


def _enrich_spacers_from_gff(spacers_list, gff_path, genome_name):
    """Parse MinCED GFF to add array coordinates and repeat info to spacer records."""
    array_info = {}
    try:
        with open(gff_path, "r") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 9 and "CRISPR" in parts[2]:
                    aid = None
                    rpt = None
                    for attr in parts[8].split(";"):
                        if attr.startswith("ID="):
                            aid = attr.split("=", 1)[1]
                        if "rpt_unit_seq=" in attr:
                            rpt = attr.split("=", 1)[1]
                    if aid:
                        array_info[aid] = {
                            "array_start": int(parts[3]),
                            "array_end": int(parts[4]),
                        }
                        if rpt:
                            array_info[aid]["repeat_sequence"] = rpt
                            array_info[aid]["repeat_length"] = len(rpt)
    except Exception:
        pass

    for spacer in spacers_list:
        if spacer["host_genome"] == genome_name:
            info = array_info.get(spacer.get("array_id"), {})
            spacer.setdefault("array_start", info.get("array_start"))
            spacer.setdefault("array_end", info.get("array_end"))
            spacer.setdefault("repeat_sequence", info.get("repeat_sequence", ""))
            spacer.setdefault("repeat_length", info.get("repeat_length", 0))


def run_minced_on_genomes(genome_files, minced_binary, progress_callback=None):
    """Run MinCED on host genome FASTA files to extract CRISPR spacers.

    Returns (spacers_df, spacers_fasta_text, summary_df).
    """
    all_spacers = []
    summaries = []
    fasta_lines = []
    total = len(genome_files)

    with tempfile.TemporaryDirectory() as tmpdir:
        for idx, uf in enumerate(genome_files):
            genome_name = uf.name.replace(".fasta", "").replace(".fa", "").replace(".fna", "")
            fasta_path = os.path.join(tmpdir, uf.name)
            with open(fasta_path, "wb") as f:
                f.write(uf.getvalue())

            txt_path = os.path.join(tmpdir, f"{uf.name}.txt")
            gff_path = os.path.join(tmpdir, f"{uf.name}.gff")
            # MinCED auto-creates spacers as <txt_basename>_spacers.fa
            spacers_path = os.path.join(tmpdir, f"{uf.name}_spacers.fa")

            # MinCED syntax: minced [options] sequence.fa [outputFile] [outputGFF]
            cmd = [minced_binary, "-spacers", "-gffFull", fasta_path, txt_path, gff_path]

            try:
                subprocess.run(cmd, capture_output=True, timeout=600)

                if os.path.isfile(spacers_path):
                    spacer_index = 0
                    for record in SeqIO.parse(spacers_path, "fasta"):
                        header_parts = record.id.split("_")
                        array_id = "_".join(header_parts[:-1]) if len(header_parts) > 1 else record.id
                        seq_str = str(record.seq)

                        all_spacers.append({
                            "host_genome": genome_name,
                            "array_id": array_id,
                            "spacer_index": spacer_index,
                            "spacer_id": f"{genome_name}__spacer_{spacer_index}",
                            "spacer_sequence": seq_str,
                            "spacer_length": len(seq_str),
                        })

                        fasta_lines.append(f">{genome_name}__spacer_{spacer_index}")
                        fasta_lines.append(seq_str)
                        spacer_index += 1

                    if os.path.isfile(gff_path):
                        _enrich_spacers_from_gff(all_spacers, gff_path, genome_name)

                    genome_spacers = [s for s in all_spacers if s["host_genome"] == genome_name]
                    unique_arrays = len(set(s["array_id"] for s in genome_spacers))
                    avg_len = np.mean([s["spacer_length"] for s in genome_spacers]) if genome_spacers else 0

                    summaries.append({
                        "host_genome": genome_name,
                        "total_arrays": unique_arrays,
                        "total_spacers": len(genome_spacers),
                        "avg_spacer_length": round(avg_len, 1),
                    })
                else:
                    summaries.append({
                        "host_genome": genome_name,
                        "total_arrays": 0, "total_spacers": 0, "avg_spacer_length": 0,
                    })
            except Exception:
                summaries.append({
                    "host_genome": genome_name,
                    "total_arrays": 0, "total_spacers": 0, "avg_spacer_length": 0,
                })

            if progress_callback:
                progress_callback((idx + 1) / total)

    spacers_df = pd.DataFrame(all_spacers) if all_spacers else pd.DataFrame()
    summary_df = pd.DataFrame(summaries) if summaries else pd.DataFrame()
    spacers_fasta = "\n".join(fasta_lines)
    return spacers_df, spacers_fasta, summary_df


def build_blast_db(source, makeblastdb_binary, db_title="plasmid_db"):
    """Build a BLAST nucleotide database.

    Args:
        source: str (file path) or list of UploadedFile objects.
        makeblastdb_binary: path to makeblastdb.

    Returns (db_path, tmpdir_handle) — caller must keep tmpdir alive.
    """
    tmpdir = tempfile.TemporaryDirectory()
    try:
        if isinstance(source, str):
            fasta_path = source
        else:
            fasta_path = os.path.join(tmpdir.name, "plasmids_combined.fasta")
            with open(fasta_path, "w") as out_f:
                for uf in source:
                    content = uf.getvalue().decode("utf-8", errors="replace")
                    if not content.endswith("\n"):
                        content += "\n"
                    out_f.write(content)

        db_path = os.path.join(tmpdir.name, db_title)
        cmd = [makeblastdb_binary, "-in", fasta_path, "-dbtype", "nucl",
               "-out", db_path, "-title", db_title, "-parse_seqids"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if result.returncode != 0:
            tmpdir.cleanup()
            return None, None
        return db_path, tmpdir
    except Exception:
        tmpdir.cleanup()
        return None, None


@st.cache_resource(show_spinner="Building reference BLAST database (one-time)...")
def get_or_build_reference_blastdb(makeblastdb_binary):
    """Build or locate pre-built BLAST database from reference plasmid sequences."""
    if os.path.isfile(REFERENCE_BLASTDB_PATH + ".ndb") or os.path.isfile(REFERENCE_BLASTDB_PATH + ".nsq"):
        return REFERENCE_BLASTDB_PATH
    if not os.path.isfile(REFERENCE_FASTA_PATH):
        return None
    os.makedirs(os.path.dirname(REFERENCE_BLASTDB_PATH), exist_ok=True)
    cmd = [makeblastdb_binary, "-in", REFERENCE_FASTA_PATH, "-dbtype", "nucl",
           "-out", REFERENCE_BLASTDB_PATH, "-title", "pLIN_reference_plasmids",
           "-parse_seqids"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode == 0:
            return REFERENCE_BLASTDB_PATH
    except Exception:
        pass
    return None


def run_spacer_blast(spacers_fasta, db_path, blastn_binary):
    """Run BLASTN-short: spacers vs plasmid BLAST database."""
    if not spacers_fasta.strip():
        return pd.DataFrame()

    with tempfile.TemporaryDirectory() as tmpdir:
        query_path = os.path.join(tmpdir, "spacers_query.fasta")
        out_path = os.path.join(tmpdir, "blast_results.tsv")

        with open(query_path, "w") as f:
            f.write(spacers_fasta)

        cmd = [
            blastn_binary, "-task", "blastn-short",
            "-query", query_path, "-db", db_path, "-out", out_path,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
            "-evalue", str(CRISPR_BLAST_EVALUE),
            "-num_threads", "4", "-max_target_seqs", "500",
            "-dust", "no", "-word_size", "7",
        ]

        try:
            subprocess.run(cmd, capture_output=True, timeout=1800)
            if os.path.isfile(out_path) and os.path.getsize(out_path) > 0:
                col_names = ["qseqid", "sseqid", "pident", "length", "mismatch",
                             "gapopen", "qstart", "qend", "sstart", "send",
                             "evalue", "bitscore", "qlen", "slen"]
                return pd.read_csv(out_path, sep="\t", header=None, names=col_names)
        except Exception:
            pass
    return pd.DataFrame()


def filter_blast_hits(blast_df):
    """Apply stringent CRISPR spacer-plasmid matching criteria."""
    if blast_df.empty:
        return pd.DataFrame()

    filtered = blast_df[
        (blast_df["gapopen"] == CRISPR_BLAST_GAPS_MAX) &
        (blast_df["pident"] >= CRISPR_BLAST_IDENTITY_MIN) &
        (blast_df["length"] >= CRISPR_BLAST_ALIGNMENT_MIN) &
        (blast_df["mismatch"] <= CRISPR_BLAST_MISMATCH_MAX)
    ].copy()

    if filtered.empty:
        return pd.DataFrame()

    filtered["host_genome"] = filtered["qseqid"].str.rsplit("__spacer_", n=1).str[0]
    filtered["plasmid_id"] = filtered["sseqid"]
    filtered["alignment_quality"] = filtered["bitscore"] / filtered["length"]
    return filtered.reset_index(drop=True)


def compute_host_probabilities(filtered_hits_df, spacer_summary_df, temperature=None):
    """Compute host probability rankings per plasmid using softmax transformation.

    Returns (host_probs_df, summary_df).
    """
    if temperature is None:
        temperature = CRISPR_SOFTMAX_TEMPERATURE
    if filtered_hits_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    pair_stats = filtered_hits_df.groupby(["host_genome", "plasmid_id"]).agg(
        spacer_hits=("qseqid", "count"),
        unique_spacers=("qseqid", "nunique"),
        avg_quality=("alignment_quality", "mean"),
        avg_identity=("pident", "mean"),
        avg_length=("length", "mean"),
        best_bitscore=("bitscore", "max"),
    ).reset_index()

    spacer_totals = spacer_summary_df[["host_genome", "total_spacers"]].copy()
    pair_stats = pair_stats.merge(spacer_totals, on="host_genome", how="left")
    pair_stats["total_spacers"] = pair_stats["total_spacers"].fillna(1)
    pair_stats["normalized_score"] = pair_stats["unique_spacers"] / pair_stats["total_spacers"].clip(lower=1)

    all_probs = []
    summaries = []

    for plasmid_id, group in pair_stats.groupby("plasmid_id"):
        scores = group["normalized_score"].values
        # Numerically stable softmax
        shifted = scores - scores.max()
        exp_scores = np.exp(shifted / temperature)
        probabilities = exp_scores / exp_scores.sum()

        group = group.copy()
        group["probability"] = probabilities
        group["rank"] = group["probability"].rank(ascending=False, method="min").astype(int)
        group = group.sort_values("rank")
        all_probs.append(group)

        top = group.iloc[0]
        confidence = "High" if top["probability"] >= 0.7 else "Medium" if top["probability"] >= 0.4 else "Low"
        summaries.append({
            "plasmid_id": plasmid_id,
            "predicted_host": top["host_genome"],
            "probability": round(float(top["probability"]), 4),
            "spacer_hits": int(top["unique_spacers"]),
            "avg_identity": round(float(top["avg_identity"]), 1),
            "confidence_category": confidence,
            "n_candidate_hosts": len(group),
        })

    host_probs_df = pd.concat(all_probs, ignore_index=True) if all_probs else pd.DataFrame()
    summary_df = pd.DataFrame(summaries) if summaries else pd.DataFrame()
    return host_probs_df, summary_df


# ══════════════════════════════════════════════════════════════════════════════
#  BACTERIAL BUDDY — OLLAMA LLM CHATBOT
# ══════════════════════════════════════════════════════════════════════════════

OLLAMA_DEFAULT_URL = "http://localhost:11434"
OLLAMA_MODELS = ["llama3.2", "llama3.1", "llama3", "mistral", "mixtral", "gemma2", "phi3"]

BACTERIAL_BUDDY_SYSTEM_PROMPT = """You are DRAGNOME Buddy, a friendly and knowledgeable AI assistant specialized in plasmid biology and antimicrobial resistance. You help researchers understand their plasmid analysis results from the pLIN (Plasmid Lineage Identification Number) classification tool.

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
                       label=f"L6 {sc} (n={sum(1 for s in strain_clusters if s == sc)})")
                       for sc in unique_strains]
    ax_dendro.legend(handles=legend_elements, title="pLIN L6",
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

    legend_elements = [mpatches.Patch(color=strain_cmap[sc], label=f"L6 {sc}") for sc in unique_strains]
    ax.legend(handles=legend_elements, title="pLIN L6", loc="lower left",
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
    ax_l.text(0.05, y, "L6 Clusters", fontsize=8, fontweight="bold", transform=ax_l.transAxes)
    y -= 0.04
    for sc in unique_strains:
        n = sum(1 for s in strain_clusters if s == sc)
        ax_l.add_patch(mpatches.FancyBboxPatch((0.05, y - 0.01), 0.08, 0.025, transform=ax_l.transAxes,
                       boxstyle="round,pad=0.003", facecolor=strain_cmap[sc], edgecolor="none"))
        ax_l.text(0.16, y + 0.003, f"L6 {sc} (n={n})", fontsize=7,
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
            "linkage_method_used", "nt_results", "prodigal_genes_df", "prodigal_summary_df",
            "mobsuite_df", "crispr_spacers_df", "crispr_spacer_summary_df",
            "crispr_host_probs_df", "crispr_host_summary_df", "crispr_filtered_hits_df",
            "mlst_df", "mlst_genome_plasmid_map", "transmission_pairs_df",
            "transmission_summary"]:
    if key not in st.session_state:
        st.session_state[key] = None
if "analysis_done" not in st.session_state:
    st.session_state.analysis_done = False
if "is_query_mode" not in st.session_state:
    st.session_state.is_query_mode = False
if "query_metadata" not in st.session_state:
    st.session_state.query_metadata = None

# DRAGNOME Buddy chat state
if "buddy_messages" not in st.session_state:
    st.session_state.buddy_messages = []
if "buddy_model" not in st.session_state:
    st.session_state.buddy_model = None


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN HEADER & UPLOAD
# ══════════════════════════════════════════════════════════════════════════════

# ── Logo header ──────────────────────────────────────────────────────────────
_logo_path = os.path.join(_LOGO_DIR, "pLIN_logo_transparent.png")
if os.path.exists(_logo_path):
    _logo_col, _title_col = st.columns([1, 4])
    with _logo_col:
        st.image(_logo_path, width=140)
    with _title_col:
        st.markdown(
            "<h1 style='margin-bottom:0; padding-top:18px;'>pLIN Classifier</h1>"
            "<p style='color:#3B6FA0; margin-top:0;'>"
            "Plasmid Lineage Identification Number System &mdash; Upload FASTA files to begin</p>",
            unsafe_allow_html=True,
        )
else:
    st.title("pLIN Classifier")
    st.caption("Plasmid Lineage Identification Number System — Upload FASTA files to begin")

# AMRFinderPlus detection (used in both upload and post-analysis views)
amr_binary, amr_db = detect_amrfinder()

# Prodigal detection
prodigal_binary = detect_prodigal()

# MOBsuite detection
mobsuite_binary = detect_mobsuite()

# Mash detection (ANI estimation)
mash_binary = detect_mash()

# FastANI detection (true ANI computation)
fastani_binary = detect_fastani()

# minimap2 detection (SNP sub-typing)
minimap2_binary = detect_minimap2()

# CRISPR Host Inference tool detection
minced_binary = detect_minced()
blastn_binary, makeblastdb_binary = detect_blastn()

# MLST chromosomal typing detection
mlst_binary = detect_mlst()

# Ollama detection for DRAGNOME Buddy
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
        metadata_file = st.file_uploader(
            "Upload metadata (optional)",
            type=["csv", "tsv", "txt"],
            accept_multiple_files=False,
            help="CSV/TSV with columns: plasmid_id (or filename), collection_date, location, patient_id, source, etc. "
                 "Used for temporal outbreak clustering and epidemiological context.",
            key="metadata_uploader",
        )
    with upload_col2:
        inc_type = st.selectbox(
            "Incompatibility Group",
            INC_GROUPS, index=0,
            help=f"'Auto-detect' uses a KNN classifier trained on {len(INC_GROUPS) - 2} Inc/Rep groups to identify Inc group per sequence",
        )
        linkage_method = st.selectbox(
            "Linkage Method",
            LINKAGE_METHODS, index=0,
            help="Single: traditional chaining (default). Complete: max distance, tighter clusters. Average: balanced. Weighted: WPGMA.",
        )
        use_adaptive = st.checkbox(
            "Adaptive thresholds",
            value=True,
            help="Calibrate pLIN thresholds per Inc group from training data distance distributions (recommended). Uncheck to use fixed universal thresholds.",
        )
        auto_filter_plasmids = st.checkbox(
            "Auto-detect plasmid contigs",
            value=True,
            help="Automatically identify and separate plasmid contigs from chromosomal sequences before pLIN assignment. "
                 "Uses sequence length, 4-mer distance to plasmid training data, and FASTA header analysis. "
                 "Chromosomal contigs are excluded from pLIN classification and shown separately. "
                 "Recommended when uploading whole-genome assemblies.",
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

        # MOBsuite mobility typing (optional)
        if mobsuite_binary:
            run_mobsuite = st.checkbox(
                "Run MOBsuite typing",
                value=False,
                help="Classify plasmid mobility, relaxase families (MOB), and MPF types using MOBsuite mob_typer.",
            )
        else:
            run_mobsuite = False
            st.info("MOBsuite not found. Install: `conda install -c bioconda mob_suite`", icon="🔬")

        # ANI Validation tools (optional)
        st.divider()
        st.markdown("**ANI Validation & SNP Sub-typing**")
        if mash_binary:
            run_mash = st.checkbox(
                "Run Mash (ANI estimation)",
                value=False,
                help="Fast MinHash-based ANI estimation. Validates pLIN cosine distances against approximate ANI values.",
            )
        else:
            run_mash = False
            st.info("Mash not found. Install: `conda install -c bioconda mash`", icon="📐")

        if fastani_binary:
            run_fastani = st.checkbox(
                "Run FastANI (true ANI)",
                value=False,
                help="Compute true Average Nucleotide Identity for all plasmid pairs. More accurate but slower than Mash.",
            )
        else:
            run_fastani = False
            st.info("FastANI not found. Install: `conda install -c bioconda fastani`", icon="📐")

        if minimap2_binary:
            run_snp_subtype = st.checkbox(
                "Run SNP sub-typing (L6 clusters)",
                value=False,
                help="Align plasmids within L6 clusters using minimap2 to count SNP differences. Critical for outbreak-level resolution.",
            )
        else:
            run_snp_subtype = False
            st.info("minimap2 not found. Install: `conda install -c bioconda minimap2`", icon="🔬")

        # CRISPR Host Inference (optional)
        st.divider()
        st.markdown("**CRISPR Host Inference**")
        crispr_tools_ok = bool(minced_binary and blastn_binary and makeblastdb_binary)
        if crispr_tools_ok:
            run_crispr = st.checkbox(
                "Run CRISPR host inference",
                value=False,
                help="Infer plasmid-host relationships using CRISPR spacer matching. "
                     "Requires host bacterial genome FASTAs.",
            )
            if run_crispr:
                crispr_source = st.radio(
                    "Plasmid database source",
                    ["Uploaded plasmids", "Reference DB (72,959 plasmids)"],
                    index=0,
                    help="Match spacers against your uploaded plasmids or the built-in reference database.",
                    horizontal=True,
                )
                crispr_host_files = st.file_uploader(
                    "Upload host bacterial genome FASTAs",
                    type=["fasta", "fa", "fna"],
                    accept_multiple_files=True,
                    help="Upload one or more bacterial genome FASTA files to screen for CRISPR spacers.",
                    key="crispr_host_uploader",
                )
            else:
                crispr_source = "Uploaded plasmids"
                crispr_host_files = None
        else:
            run_crispr = False
            crispr_source = "Uploaded plasmids"
            crispr_host_files = None
            missing = []
            if not minced_binary:
                missing.append("minced")
            if not blastn_binary or not makeblastdb_binary:
                missing.append("blast")
            st.info(
                f"CRISPR Host Inference requires: {', '.join(missing)}. "
                f"Install: `conda install -c bioconda {' '.join(missing)}`",
                icon="🧫",
            )

        # ── Chromosomal Typing (MLST) ────────────────────────────────
        st.divider()
        st.markdown("**Chromosomal Typing (MLST)**")
        if mlst_binary:
            run_mlst = st.checkbox(
                "Run MLST typing",
                value=False,
                help="Type host bacterial genomes using MLST (Torsten Seemann's mlst tool). "
                     "Combines with pLIN to distinguish clonal vs horizontal plasmid spread. "
                     "Uses the same host genome FASTAs as CRISPR analysis if available.",
            )
            if run_mlst:
                if crispr_tools_ok and run_crispr and crispr_host_files:
                    mlst_genome_files = crispr_host_files
                    st.caption(f"Using {len(mlst_genome_files)} genome(s) from CRISPR uploader")
                else:
                    mlst_genome_files = st.file_uploader(
                        "Upload host bacterial genome FASTAs",
                        type=["fasta", "fa", "fna"],
                        accept_multiple_files=True,
                        help="Upload assembled bacterial genome FASTAs for MLST typing.",
                        key="mlst_genome_uploader",
                    )
            else:
                mlst_genome_files = None
        else:
            run_mlst = False
            mlst_genome_files = None
            st.info(
                "mlst not found. Install: `conda install -c bioconda mlst`",
                icon="🧬",
            )

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
    auto_filter_plasmids = st.session_state.get("_auto_filter_plasmids", True)
    use_nt = st.session_state.get("_use_nt", False)
    nt_model_choice = st.session_state.get("_nt_model_choice", None)
    nt_device = st.session_state.get("_nt_device", "cpu")
    run_amr = amr_binary is not None
    run_prodigal = st.session_state.get("_run_prodigal", False)
    run_mash = False
    run_fastani = False
    run_snp_subtype = False
    metadata_file = None
    run_btn = False


# ── Sidebar ───────────────────────────────────────────────────────────────────

with st.sidebar:
    if os.path.exists(_logo_path):
        st.image(_logo_path, width=120)
    else:
        st.title("pLIN")
    st.caption("Plasmid Lineage Identification Number")
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
    st.session_state._auto_filter_plasmids = auto_filter_plasmids
    st.session_state._use_nt = use_nt
    st.session_state._nt_model_choice = nt_model_choice
    st.session_state._nt_device = nt_device
    st.session_state._run_prodigal = run_prodigal

    # Parse metadata CSV if provided
    metadata_df = None
    if metadata_file is not None:
        try:
            sep = "\t" if metadata_file.name.endswith((".tsv", ".txt")) else ","
            metadata_df = pd.read_csv(metadata_file, sep=sep)
            # Standardize join column: try plasmid_id, filename, sample_id
            join_col = None
            for candidate in ["plasmid_id", "Plasmid_ID", "filename", "Filename", "sample_id", "Sample_ID", "name", "Name"]:
                if candidate in metadata_df.columns:
                    join_col = candidate
                    break
            if join_col and join_col != "plasmid_id":
                metadata_df = metadata_df.rename(columns={join_col: "plasmid_id"})
            # Parse date columns
            for col in metadata_df.columns:
                if any(kw in col.lower() for kw in ["date", "collection_date", "sample_date"]):
                    metadata_df[col] = pd.to_datetime(metadata_df[col], errors="coerce")
            st.session_state.metadata_df = metadata_df
        except Exception as e:
            st.warning(f"Could not parse metadata file: {e}", icon="⚠️")
            metadata_df = None
    else:
        st.session_state.metadata_df = None

    progress = st.progress(0, text="Starting analysis...")

    # Step 1: Parse sequences (+ Inc group auto-detection)
    progress.progress(5, text="Parsing FASTA files & detecting Inc groups..." if inc_type == "Auto-detect"
                      else "Parsing FASTA files...")
    try:
        records = parse_uploaded_fastas(uploaded_files, inc_type)
        if len(records) < 1:
            st.error("No valid sequences found in uploaded files.")
            st.stop()
        st.session_state.records = records
        # Run duplicate detection on loaded records
        if len(records) >= 2:
            st.session_state.duplicate_pairs = detect_duplicates(records)
        else:
            st.session_state.duplicate_pairs = []
    except Exception as e:
        st.error(f"Failed to parse FASTA files: {e}")
        st.stop()

    # Step 2: K-mer vectors
    progress.progress(15, text=f"Computing 4-mer vectors for {len(records)} sequences...")
    sequences = [r["sequence"] for r in records]
    vectors = compute_kmer_vectors(tuple(sequences), k=4)

    # Validate vectors shape matches records count
    if vectors.shape[0] != len(records):
        st.error(
            f"Vector computation mismatch: {len(records)} records but {vectors.shape[0]} vectors. "
            "This is likely a caching issue. Please refresh the page (Ctrl+R / Cmd+R) and try again."
        )
        st.stop()

    # Step 2a: Automatic plasmid vs chromosome classification
    original_records = list(records)  # preserve pre-filter list for strain summary
    if auto_filter_plasmids:
        progress.progress(20, text="Classifying contigs: plasmid vs chromosome...")
        contig_classes = classify_contigs_plasmid_vs_chromosome(records, vectors)
        st.session_state.contig_classification = contig_classes

        # Three-way split: plasmid (>=95% conf), incomplete_plasmid, chromosome
        plasmid_indices = [i for i, c in enumerate(contig_classes) if c["classification"] == "plasmid"]
        incomplete_indices = [i for i, c in enumerate(contig_classes) if c["classification"] == "incomplete_plasmid"]
        chromo_indices = [i for i, c in enumerate(contig_classes) if c["classification"] == "chromosome"]

        # Store excluded chromosomes
        if chromo_indices:
            chromo_ids = [contig_classes[i]["plasmid_id"] for i in chromo_indices]
            chromo_df = pd.DataFrame([contig_classes[i] for i in chromo_indices])
            st.session_state.excluded_chromosomes = chromo_df
            st.info(
                f"**Excluded {len(chromo_indices)} chromosomal contig(s)**: "
                f"{', '.join(chromo_ids[:5])}{'...' if len(chromo_ids) > 5 else ''}. "
                f"These are too distant from known plasmid groups or too large to be plasmids.",
                icon="🧬",
            )
        else:
            st.session_state.excluded_chromosomes = None

        # Store incomplete plasmids — these get AMR/mobility/other modules
        # but NOT pLIN code assignment
        if incomplete_indices:
            incomplete_ids = [contig_classes[i]["plasmid_id"] for i in incomplete_indices]
            incomplete_df = pd.DataFrame([contig_classes[i] for i in incomplete_indices])
            st.session_state.incomplete_plasmids = incomplete_df
            # Keep incomplete plasmids in records for AMR/mobility analysis
            # but mark them so pLIN assignment is skipped
            for idx in incomplete_indices:
                records[idx]["_skip_plin"] = True
            st.warning(
                f"**{len(incomplete_indices)} contig(s) classified as incomplete/uncertain plasmid** "
                f"(confidence <95%): {', '.join(incomplete_ids[:5])}"
                f"{'...' if len(incomplete_ids) > 5 else ''}. "
                f"These will be analysed for AMR genes, mobility, and other features "
                f"but will **not** receive pLIN code assignment.",
                icon="⚠️",
            )
        else:
            st.session_state.incomplete_plasmids = None

        # Only truly chromosomal contigs are removed entirely
        non_chromo_indices = plasmid_indices + incomplete_indices
        non_chromo_indices.sort()

        if not non_chromo_indices:
            st.error(
                "No plasmid contigs detected. All uploaded sequences appear to be chromosomal. "
                "If this is incorrect, uncheck 'Auto-detect plasmid contigs' in the sidebar "
                "to force pLIN assignment on all sequences."
            )
            st.stop()

        # Keep both plasmid and incomplete_plasmid contigs for downstream analysis
        records = [records[i] for i in non_chromo_indices]
        vectors = vectors[non_chromo_indices]
        st.session_state.records = records

        # Track which contigs should get pLIN codes (pre-merge IDs)
        plin_eligible_original = set(
            contig_classes[i]["plasmid_id"] for i in plasmid_indices
        )

        # ── Strain-aware assembly summary ─────────────────────────────
        # Build per-source-file contig summary for multi-contig uploads
        source_files = {}
        for idx, cc in enumerate(contig_classes):
            src = original_records[idx].get("_strain_id", original_records[idx].get("source_file", "unknown"))
            if src not in source_files:
                source_files[src] = {"plasmid": [], "incomplete_plasmid": [], "chromosome": []}
            source_files[src][cc["classification"]].append(cc)

        # Only show assembly summary when any file has >1 contig
        multi_contig_files = {sf: cats for sf, cats in source_files.items()
                              if sum(len(v) for v in cats.values()) > 1}
        if multi_contig_files:
            summary_rows = []
            for sf, cats in multi_contig_files.items():
                n_plas = len(cats["plasmid"])
                n_inc = len(cats["incomplete_plasmid"])
                n_chr = len(cats["chromosome"])
                n_total = n_plas + n_inc + n_chr
                plas_ids = [c["plasmid_id"] for c in cats["plasmid"]]
                chr_sizes = [f"{c['length_bp']/1e6:.1f} Mb" for c in cats["chromosome"]]
                summary_rows.append({
                    "Source File": sf,
                    "Total Contigs": n_total,
                    "Plasmid": n_plas,
                    "Incomplete": n_inc,
                    "Chromosome": n_chr,
                    "Plasmid Contigs": ", ".join(plas_ids[:5]) + ("..." if len(plas_ids) > 5 else ""),
                })
                # Show per-file breakdown
                parts = []
                if n_plas:
                    parts.append(f"**{n_plas}** plasmid contig(s)")
                if n_inc:
                    parts.append(f"**{n_inc}** incomplete")
                if n_chr:
                    parts.append(f"**{n_chr}** chromosome ({', '.join(chr_sizes)})")
                st.info(
                    f"**Assembly** `{sf}` ({n_total} contigs): {' | '.join(parts)}",
                    icon="📋",
                )

            st.session_state.strain_contig_summary = pd.DataFrame(summary_rows)
        else:
            st.session_state.strain_contig_summary = None

        # ── Step 2c: Merge multi-contig plasmids from same file + Inc type ──
        pre_merge_count = len(records)
        records, merge_map = merge_multicontig_plasmids(records)
        st.session_state.merge_map = merge_map

        if len(records) < pre_merge_count:
            # Merging happened — recompute 4-mer vectors on merged sequences
            sequences = tuple(r["sequence"] for r in records)
            vectors = compute_kmer_vectors(sequences, k=4)

            # Show merge summary
            merged_groups = {k: v for k, v in merge_map.items() if len(v) > 1}
            for mid, orig_ids in merged_groups.items():
                st.info(
                    f"**Merged {len(orig_ids)} contigs** from same assembly + Inc type "
                    f"into single plasmid: {', '.join(orig_ids)}",
                    icon="🔗",
                )

        st.session_state.records = records

        # Update pLIN eligible IDs to use merged record IDs
        plin_eligible = set()
        for rec in records:
            merged_from = rec.get("_merged_from")
            if merged_from:
                # Merged record is eligible if ANY original contig was eligible
                if any(oid in plin_eligible_original for oid in merged_from):
                    plin_eligible.add(rec["plasmid_id"])
            else:
                if rec["plasmid_id"] in plin_eligible_original:
                    plin_eligible.add(rec["plasmid_id"])
        st.session_state._plin_eligible_ids = plin_eligible

    else:
        st.session_state.contig_classification = None
        st.session_state.excluded_chromosomes = None
        st.session_state.incomplete_plasmids = None
        st.session_state._plin_eligible_ids = None
        st.session_state.strain_contig_summary = None
        st.session_state.merge_map = None

    # Determine analysis mode
    is_query_mode = len(records) == 1

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
    if is_query_mode:
        # Single-plasmid query mode: nearest-neighbour lookup against training DB
        progress.progress(40, text="Assigning pLIN code via nearest-neighbour lookup...")
        try:
            plin_codes, cluster_assignments, query_metadata = assign_plin_query_mode(
                vectors, records, thresholds=active_thresholds
            )
            Z = None
            dist_condensed = None
            st.session_state.query_metadata = query_metadata
        except FileNotFoundError as e:
            st.error(str(e))
            st.stop()
    else:
        # Multi-plasmid mode: pLIN codes are reference-anchored (each contig
        # inherits its code from its own nearest reference-database neighbour),
        # so results are reproducible and consistent with single-contig query
        # mode and with previously published/reference-based runs. De novo
        # pairwise clustering (Z/dist_condensed) is still computed separately,
        # purely to drive the dendrogram/heatmap visualisation of relatedness
        # among the uploaded contigs themselves.
        progress.progress(40, text="Assigning pLIN codes via nearest-neighbour lookup...")
        try:
            plin_codes, cluster_assignments, query_metadata = assign_plin_query_mode(
                vectors, records, thresholds=active_thresholds
            )
            st.session_state.query_metadata = query_metadata
        except FileNotFoundError as e:
            st.error(str(e))
            st.stop()

        progress.progress(45, text=f"Clustering ({linkage_method} linkage) for visualisation...")
        _, _, Z, dist_condensed = assign_plin_codes(
            vectors, linkage_method=linkage_method, thresholds=active_thresholds
        )

    strain_clusters = list(cluster_assignments["F"])
    st.session_state.linkage_method_used = linkage_method
    st.session_state.is_query_mode = is_query_mode

    # Step 4: Build results
    progress.progress(55, text="Building results table...")
    plin_df = build_results_df(records, plin_codes, cluster_assignments)
    # Use plasmid_id for labels (unique per sequence) instead of source_file
    labels = [r["plasmid_id"] for r in records]

    # Add query-mode metadata columns if available
    if st.session_state.query_metadata:
        qm = st.session_state.query_metadata
        plin_df["nn_plasmid"] = [m["nn_plasmid"] for m in qm]
        plin_df["nn_distance"] = [m["nn_distance"] for m in qm]
        plin_df["nn_inc_type"] = [m["nn_inc_type"] for m in qm]
        plin_df["nn_plin"] = [m["nn_plin"] for m in qm]

    # Mark incomplete plasmids — they keep AMR/mobility results but no pLIN code
    plin_eligible = st.session_state.get("_plin_eligible_ids")
    if plin_eligible is not None:
        for idx, row in plin_df.iterrows():
            pid = row.get("plasmid_id", row.get("Plasmid", ""))
            if pid not in plin_eligible:
                plin_df.at[idx, "pLIN"] = "N/A (incomplete plasmid)"
                plin_df.at[idx, "inc_type"] = plin_df.at[idx, "inc_type"] + " *"
                # Blank out hierarchical level assignments
                for lvl in ["A", "B", "C", "D", "E", "F"]:
                    if lvl in plin_df.columns:
                        plin_df.at[idx, lvl] = "—"

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

    # Step 5c: MOBsuite mobility typing (optional)
    mobsuite_df = pd.DataFrame()
    if run_mobsuite and mobsuite_binary:
        progress.progress(86, text="Running MOBsuite mobility typing...")

        def mob_cb(pct):
            progress.progress(int(86 + pct * 4), text=f"MOBsuite: {int(pct * 100)}%")

        mobsuite_df = run_mobsuite_on_files(uploaded_files, mobsuite_binary, mob_cb)

    st.session_state.mobsuite_df = mobsuite_df

    # Step 5d: CRISPR host inference (optional)
    if run_crispr and crispr_host_files:
        progress.progress(87, text="Running CRISPR host inference — extracting spacers...")

        # 5d-i: Run MinCED on host genomes
        def crispr_cb(pct):
            progress.progress(int(87 + pct * 2), text=f"MinCED spacer extraction: {int(pct * 100)}%")

        spacers_df, spacers_fasta_text, spacer_summary_df = run_minced_on_genomes(
            crispr_host_files, minced_binary, crispr_cb
        )
        st.session_state.crispr_spacers_df = spacers_df
        st.session_state.crispr_spacer_summary_df = spacer_summary_df

        if len(spacers_df) > 0 and spacers_fasta_text.strip():
            # 5d-ii: Build or load BLAST DB
            progress.progress(90, text="Building BLAST database for spacer matching...")
            tmpdir_handle = None
            if crispr_source == "Reference DB (72,959 plasmids)":
                db_path = get_or_build_reference_blastdb(makeblastdb_binary)
            else:
                db_path, tmpdir_handle = build_blast_db(uploaded_files, makeblastdb_binary)

            if db_path:
                # 5d-iii: Run BLASTN-short
                progress.progress(91, text="Running BLASTN-short: spacers vs plasmids...")
                blast_df = run_spacer_blast(spacers_fasta_text, db_path, blastn_binary)

                if len(blast_df) > 0:
                    # 5d-iv: Filter hits
                    progress.progress(92, text="Filtering BLAST hits (stringent criteria)...")
                    filtered_hits_df = filter_blast_hits(blast_df)
                    st.session_state.crispr_filtered_hits_df = filtered_hits_df

                    if len(filtered_hits_df) > 0:
                        # 5d-v: Compute host probabilities
                        progress.progress(93, text="Computing host-plasmid probability rankings...")
                        host_probs_df, host_summary_df = compute_host_probabilities(
                            filtered_hits_df, spacer_summary_df
                        )
                        st.session_state.crispr_host_probs_df = host_probs_df
                        st.session_state.crispr_host_summary_df = host_summary_df
                    else:
                        st.session_state.crispr_filtered_hits_df = pd.DataFrame()
                else:
                    st.session_state.crispr_filtered_hits_df = pd.DataFrame()

                # Clean up temp BLAST DB if built from uploads
                if tmpdir_handle is not None:
                    try:
                        tmpdir_handle.cleanup()
                    except Exception:
                        pass
            else:
                st.warning("Could not build BLAST database for CRISPR analysis.")
        else:
            st.info("No CRISPR spacers were extracted from the uploaded host genomes.")

    # Step 5e: MLST chromosomal typing (optional)
    if run_mlst and mlst_genome_files:
        progress.progress(93, text="Running MLST on host genomes...")

        def mlst_cb(pct):
            progress.progress(int(93 + pct * 1), text=f"MLST typing: {int(pct * 100)}%")

        mlst_df = run_mlst_on_genomes(mlst_genome_files, mlst_binary, mlst_cb)
        st.session_state.mlst_df = mlst_df

        if len(mlst_df) > 0:
            plasmid_ids = plin_df["plasmid_id"].tolist()
            genome_plasmid_map = build_genome_plasmid_mapping(
                mlst_genome_files, plasmid_ids,
                metadata_df=st.session_state.get("metadata_df"),
            )
            st.session_state.mlst_genome_plasmid_map = genome_plasmid_map

            if len(genome_plasmid_map) > 0:
                transmission_pairs_df, transmission_summary = classify_transmission_mode(
                    mlst_df, plin_df, genome_plasmid_map,
                )
                st.session_state.transmission_pairs_df = transmission_pairs_df
                st.session_state.transmission_summary = transmission_summary

    # Step 6: Integration
    progress.progress(94, text="Integrating results...")
    integrated_df = integrate_plin_amr(plin_df, amr_df)

    # Merge user-provided metadata if available
    if st.session_state.get("metadata_df") is not None:
        meta = st.session_state.metadata_df
        if "plasmid_id" in meta.columns:
            # Strip file extensions from metadata plasmid_id to match
            meta_clean = meta.copy()
            meta_clean["plasmid_id"] = meta_clean["plasmid_id"].astype(str).str.replace(
                r"\.(fasta|fa|fna)$", "", regex=True
            )
            # Merge, keeping all plasmids (left join)
            extra_cols = [c for c in meta_clean.columns if c != "plasmid_id" and c not in integrated_df.columns]
            if extra_cols:
                integrated_df = integrated_df.merge(
                    meta_clean[["plasmid_id"] + extra_cols], on="plasmid_id", how="left"
                )

    st.session_state.integrated_df = integrated_df

    # Step 7: Mobility prediction
    progress.progress(95, text="Predicting plasmid mobility...")
    mobsuite_df = st.session_state.get("mobsuite_df")
    mobility_results = []
    for rec in records:
        mob_result = classify_mobility(amr_df, rec["source_file"], mobsuite_df=mobsuite_df)
        mobility_results.append({
            "plasmid_id": rec["plasmid_id"],
            "mobility": mob_result["mobility"],
            "mobility_genes": "; ".join(mob_result["genes"]) if mob_result["genes"] else "",
            "mobility_detail": mob_result["detail"],
            "mobility_source": mob_result["source"],
            "relaxase_family": mob_result["relaxase_family"],
            "mpf_type": mob_result["mpf_type"],
        })
    st.session_state.mobility_results = pd.DataFrame(mobility_results)

    # Step 8: Outbreak detection
    progress.progress(92, text="Scanning for outbreak clusters...")
    outbreak_clusters = detect_outbreak_clusters(plin_df, integrated_df)
    st.session_state.outbreak_clusters = outbreak_clusters

    # Step 8b: Temporal outbreak clustering (if metadata with dates provided)
    temporal_clusters = []
    if st.session_state.get("metadata_df") is not None:
        progress.progress(93, text="Checking temporal outbreak patterns...")
        temporal_clusters = detect_temporal_outbreak_clusters(
            plin_df, integrated_df, st.session_state.metadata_df, time_window_days=30
        )
    st.session_state.temporal_outbreak_clusters = temporal_clusters

    # Step 9: Mash ANI estimation (optional)
    if run_mash and mash_binary:
        progress.progress(94, text="Running Mash ANI estimation...")
        mash_df = run_mash_distances(uploaded_files, mash_binary)
        st.session_state.mash_df = mash_df
        # ANI concordance check — cross-validate pLIN cosine distances vs Mash ANI
        st.session_state.concordance_df = check_ani_concordance(plin_df, mash_df)
    else:
        st.session_state.mash_df = None
        st.session_state.concordance_df = pd.DataFrame()

    # Step 10: FastANI (optional)
    if run_fastani and fastani_binary:
        progress.progress(96, text="Running FastANI (true ANI computation)...")
        fastani_df = run_fastani(uploaded_files, fastani_binary)
        st.session_state.fastani_df = fastani_df
    else:
        st.session_state.fastani_df = None

    # Step 11: SNP sub-typing within L6 clusters (optional)
    if run_snp_subtype and minimap2_binary:
        progress.progress(97, text="SNP sub-typing within L6 clusters...")
        snp_df = run_snp_subtyping(records, cluster_assignments, minimap2_binary)
        st.session_state.snp_subtype_df = snp_df

        # L5: Recombination detection (runs alongside SNP subtyping)
        progress.progress(98, text="Detecting recombination signals...")
        recom_df = detect_recombination_signals(records, cluster_assignments, minimap2_binary)
        st.session_state.recombination_df = recom_df

        # L7: Evolutionary rate estimation (needs SNP data + dates)
        metadata_df = st.session_state.get("metadata_df")
        if metadata_df is not None and snp_df is not None:
            evo_df = estimate_evolutionary_rate(snp_df, metadata_df, records)
            st.session_state.evo_rate_df = evo_df
        else:
            st.session_state.evo_rate_df = pd.DataFrame()
    else:
        st.session_state.snp_subtype_df = None
        st.session_state.recombination_df = pd.DataFrame()
        st.session_state.evo_rate_df = pd.DataFrame()

    # L8: Cluster stability assessment (if >=4 plasmids and de novo mode)
    if len(records) >= 4 and st.session_state.Z is not None:
        progress.progress(98, text="Assessing cluster stability...")
        stability_df = assess_cluster_stability(vectors, plin_codes, n_bootstrap=50)
        st.session_state.stability_df = stability_df
        linkage_df = compare_linkage_methods(vectors)
        st.session_state.linkage_comparison_df = linkage_df
    else:
        st.session_state.stability_df = pd.DataFrame()
        st.session_state.linkage_comparison_df = pd.DataFrame()

    # ── L3: Assembly completeness assessment ────────────────────────────────
    progress.progress(99, text="Assessing assembly completeness...")
    completeness_df = assess_assembly_completeness(
        records, prodigal_summary_df=prodigal_summary_df
    )
    st.session_state.completeness_df = completeness_df

    # ── L4: Database coverage assessment ──────────────────────────────────
    nn_dists = None
    inc_types_list = plin_df["predicted_inc"].tolist() if "predicted_inc" in plin_df.columns else []
    if "nn_distance" in plin_df.columns:
        nn_dists = plin_df["nn_distance"].tolist()
    if len(inc_types_list) > 0:
        coverage_df = assess_database_coverage(vectors, inc_types_list, nn_dists)
        st.session_state.coverage_df = coverage_df
    else:
        st.session_state.coverage_df = pd.DataFrame()

    # ── L6: Novel Inc group discovery ─────────────────────────────────────
    if "inc_confidence" in plin_df.columns:
        conf_scores = plin_df["inc_confidence"].tolist()
        novel_df = discover_novel_inc_groups(vectors, inc_types_list, conf_scores, records)
        st.session_state.novel_inc_df = novel_df
    else:
        st.session_state.novel_inc_df = pd.DataFrame()

    # ── L10: MGE boundary detection (if Prodigal results available) ───────
    if len(prodigal_genes_df) > 0:
        mge_df = detect_mge_boundaries(prodigal_genes_df)
        st.session_state.mge_df = mge_df
    else:
        st.session_state.mge_df = pd.DataFrame()

    progress.progress(100, text="Analysis complete!")
    st.session_state.analysis_done = True
    st.rerun()


# ══════════════════════════════════════════════════════════════════════════════
#  TABS
# ══════════════════════════════════════════════════════════════════════════════

tab_overview, tab_results, tab_clado, tab_amr, tab_epi, tab_crispr, tab_buddy, tab_export = st.tabs(
    ["📋 Overview", "📊 Results", "🌳 Cladogram", "💊 AMR Analysis",
     "🔬 Epidemiology", "🧫 CRISPR Host", "🧬 DRAGNOME Buddy", "📥 Export"]
)

# ── TAB 1: Overview ──────────────────────────────────────────────────────────

with tab_overview:
    st.header("pLIN — Plasmid Lineage Identification Number")
    st.markdown("""
    **pLIN** assigns each plasmid a six-position hierarchical code (`L1.L2.L3.L4.L5.L6`)
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
        2. **Auto-detect** Inc/Rep group (KNN classifier, 91.1% accuracy, 28 groups)
        3. **Compute** tetranucleotide frequency vectors (256 features)
        4. **Calculate** pairwise cosine distances
        5. **Cluster** using single-linkage hierarchical method
        6. **Cut** tree at 6 thresholds → pLIN codes
        7. **Screen** for AMR genes with AMRFinderPlus *(optional)*
        """)

    if st.session_state.analysis_done:
        if st.session_state.get("is_query_mode"):
            st.success(f"Query mode — pLIN assigned via nearest-neighbour lookup against {len(st.session_state.get('X_ref', []))}-plasmid training database.")
        else:
            st.success("Analysis complete!")
        df = st.session_state.plin_df
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Plasmids", len(df))
        c2.metric("Unique pLIN Codes", df["pLIN"].nunique())
        c3.metric("L6 Clusters (F)", df["bin_F"].nunique())
        amr_count = st.session_state.integrated_df["AMR_count"].sum() if st.session_state.integrated_df is not None else 0
        c4.metric("AMR Detections", int(amr_count))

        # ── Input Sequence Quality Report ─────────────────────────────────
        records = st.session_state.get("records", [])
        all_quality_warnings = []
        for rec in records:
            for w in rec.get("quality_warnings", []):
                all_quality_warnings.append(w)

        n_errors = sum(1 for w in all_quality_warnings if w["level"] == "error")
        n_warnings = sum(1 for w in all_quality_warnings if w["level"] == "warning")
        n_clean = len(records) - len({w["message"].split(":")[0] for w in all_quality_warnings if w["level"] in ("error", "warning")})

        if all_quality_warnings:
            st.subheader("Input Sequence Quality")
            qc1, qc2, qc3 = st.columns(3)
            qc1.metric("Passed", n_clean, delta=None)
            qc2.metric("Warnings", n_warnings, delta=None)
            qc3.metric("Errors", n_errors, delta=None)

            if n_errors > 0:
                error_msgs = [w for w in all_quality_warnings if w["level"] == "error"]
                st.error(f"**{n_errors} sequence quality error(s)** — results for these sequences are unreliable.")
                with st.expander(f"View errors ({n_errors})"):
                    st.dataframe(
                        pd.DataFrame(error_msgs)[["check", "message", "value"]],
                        use_container_width=True, hide_index=True,
                    )

            if n_warnings > 0:
                warn_msgs = [w for w in all_quality_warnings if w["level"] == "warning"]
                st.warning(f"**{n_warnings} sequence quality warning(s)** — interpret flagged assignments with caution.")
                with st.expander(f"View warnings ({n_warnings})"):
                    st.dataframe(
                        pd.DataFrame(warn_msgs)[["check", "message", "value"]],
                        use_container_width=True, hide_index=True,
                    )

            # Duplicate detection
            duplicate_pairs = st.session_state.get("duplicate_pairs", [])
            if duplicate_pairs:
                st.info(f"**{len(duplicate_pairs)} near-identical sequence pair(s) detected** — these may be the same plasmid uploaded twice.")
                with st.expander(f"View duplicate pairs ({len(duplicate_pairs)})"):
                    dup_df = pd.DataFrame(duplicate_pairs, columns=["Plasmid A", "Plasmid B", "Cosine Distance"])
                    st.dataframe(dup_df, use_container_width=True, hide_index=True)
        else:
            if len(records) > 0:
                st.success(f"All {len(records)} sequences passed input quality checks.")

        # Sequence length warning for short plasmids
        short_plasmids = df[df["length_bp"] < SHORT_PLASMID_THRESHOLD]
        if len(short_plasmids) > 0:
            st.warning(
                f"**{len(short_plasmids)} plasmid(s) shorter than {SHORT_PLASMID_THRESHOLD/1000:.0f} kb** — "
                f"4-mer frequency vectors from short sequences have higher stochastic variance, "
                f"which may produce unreliable pLIN codes and inflated inter-plasmid distances. "
                f"Interpret these assignments with caution."
            )
            with st.expander(f"View short plasmids ({len(short_plasmids)})"):
                st.dataframe(
                    short_plasmids[["plasmid_id", "length_bp", "pLIN", "inc_type"]],
                    use_container_width=True, hide_index=True,
                )

        # Show excluded chromosomal contigs (from auto-detection)
        excluded_chromo = st.session_state.get("excluded_chromosomes")
        if excluded_chromo is not None and len(excluded_chromo) > 0:
            st.success(
                f"**{len(excluded_chromo)} chromosomal contig(s) auto-detected and excluded** — "
                f"only plasmid sequences were assigned pLIN codes."
            )
            with st.expander(f"View excluded chromosomal contigs ({len(excluded_chromo)})"):
                st.dataframe(
                    excluded_chromo[["plasmid_id", "length_bp", "confidence", "reason"]],
                    use_container_width=True, hide_index=True,
                )

        # Show incomplete plasmids (included for AMR/mobility but no pLIN)
        incomplete_plasm = st.session_state.get("incomplete_plasmids")
        if incomplete_plasm is not None and len(incomplete_plasm) > 0:
            st.warning(
                f"**{len(incomplete_plasm)} incomplete/uncertain plasmid contig(s)** — "
                f"analysed for AMR genes, mobility, and other features but "
                f"**not** assigned pLIN codes (confidence <95%). "
                f"These appear in results as 'N/A (incomplete plasmid)'."
            )
            with st.expander(f"View incomplete plasmid contigs ({len(incomplete_plasm)})"):
                st.dataframe(
                    incomplete_plasm[["plasmid_id", "length_bp", "confidence", "reason"]],
                    use_container_width=True, hide_index=True,
                )

        # ── Multi-Contig Assembly Summary ─────────────────────────────
        strain_summary = st.session_state.get("strain_contig_summary")
        if strain_summary is not None and len(strain_summary) > 0:
            with st.expander(f"Multi-Contig Assembly Summary ({len(strain_summary)} assemblies)"):
                st.markdown(
                    "When uploading whole-genome assemblies containing both chromosomal and "
                    "plasmid contigs, pLIN automatically identifies and separates plasmid "
                    "sequences for classification. Only confirmed plasmid contigs receive "
                    "pLIN codes."
                )
                st.dataframe(
                    strain_summary,
                    use_container_width=True, hide_index=True,
                )

        # Chromosomal sequence warning for large contigs still in results
        # (shown when auto-filter is OFF)
        large_contigs = df[df["length_bp"] > CHROMOSOMAL_THRESHOLD]
        if len(large_contigs) > 0:
            st.warning(
                f"**{len(large_contigs)} sequence(s) larger than {CHROMOSOMAL_THRESHOLD/1000:.0f} kb in results** — "
                f"these are likely chromosomal. Enable 'Auto-detect plasmid contigs' in the sidebar "
                f"to automatically exclude chromosomal sequences before pLIN assignment."
            )
            with st.expander(f"View putative chromosomal sequences ({len(large_contigs)})"):
                st.dataframe(
                    large_contigs[["plasmid_id", "length_bp", "pLIN", "inc_type"]],
                    use_container_width=True, hide_index=True,
                )

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

            # Check for multiple Inc type detections (multi-replicon plasmids)
            has_multi_inc = "inc_is_multiple" in df.columns
            if has_multi_inc:
                multi_inc_df = df[df["inc_is_multiple"] == True]
                multi_inc_count = len(multi_inc_df)
                if multi_inc_count > 0:
                    blast_count = int(df.get("inc_blast_used", False).sum()) if "inc_blast_used" in df.columns else 0
                    blast_note = f" ({blast_count} confirmed by PlasmidFinder BLAST)" if blast_count > 0 else ""
                    st.warning(
                        f"**{multi_inc_count} plasmid(s) detected with Multiple Inc types{blast_note}** — "
                        f"These are multi-replicon plasmids carrying replicons from multiple incompatibility groups. "
                        f"KNN threshold: ≥{MULTI_INC_THRESHOLD*100:.0f}% for secondary Inc types. "
                        f"For large plasmids (>100 kb) with 100% KNN confidence, PlasmidFinder BLAST is used as a secondary signal."
                    )
                    with st.expander(f"View Multiple Inc type plasmids ({multi_inc_count})"):
                        multi_cols = ["plasmid_id", "inc_type", "inc_confidence"]
                        if "inc_multiple_types" in df.columns:
                            multi_cols.append("inc_multiple_types")
                        if "inc_blast_used" in df.columns:
                            multi_cols.append("inc_blast_used")
                        if "inc_blast_hits" in df.columns:
                            multi_cols.append("inc_blast_hits")
                        multi_display = multi_inc_df[multi_cols].copy()
                        # Format blast hits for display
                        if "inc_blast_hits" in multi_display.columns:
                            multi_display["inc_blast_hits"] = multi_display["inc_blast_hits"].apply(
                                lambda hits: ", ".join(f"{g} ({p:.0f}% id)" for g, p, c in hits[:3]) if isinstance(hits, list) and hits else "—"
                            )
                        if "inc_blast_used" in multi_display.columns:
                            multi_display["inc_blast_used"] = multi_display["inc_blast_used"].apply(
                                lambda x: "✔ BLAST" if x else "KNN"
                            )
                        multi_display = multi_display.rename(columns={
                            "inc_type": "Reported Inc Type(s)",
                            "inc_confidence": "KNN Confidence",
                            "inc_multiple_types": "All Detected (KNN)",
                            "inc_blast_used": "Detection Method",
                            "inc_blast_hits": "PlasmidFinder BLAST Hits",
                        })
                        st.dataframe(multi_display, use_container_width=True, hide_index=True)
                        st.caption(
                            "**Multi-replicon plasmids** carry replicons from multiple incompatibility groups on the same molecule. "
                            "Common in large conjugative plasmids (e.g. IncHI2 + IncN co-occurring on VIM/NDM resistance plasmids). "
                            "Detection method: KNN probability split ≥15% triggers multi-replicon flag; "
                            "PlasmidFinder BLAST (enterobacteriales.fsa, ≥80% identity, ≥60% replicon coverage) "
                            "is used as a secondary signal for large plasmids where KNN gives 100% to one group. "
                            "**The pLIN code is unaffected** — it is based on whole-sequence 4-mer composition and correctly "
                            "identifies these plasmids as a single lineage regardless of which replicon label is assigned."
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
            # Don't count Unknown/Novel or Multiple Inc types as single Inc groups for the warning
            real_inc_groups = [g for g in df["inc_type"].unique()
                              if g != "Unknown/Novel" and not g.startswith("Multiple:")]
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

        # Add query-mode columns if available
        if "nn_plasmid" in df.columns:
            display_cols += [c for c in ["nn_plasmid", "nn_distance", "nn_inc_type", "nn_plin"]
                            if c in df.columns and c not in display_cols]

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

        # ── Classification Quality Report (from CV metrics) ───────────────
        cv_data = load_cv_metrics()
        if cv_data["accuracy"] > 0:
            with st.expander("Classification Quality (Cross-Validation Metrics)"):
                st.markdown(
                    f"**Overall KNN accuracy:** {cv_data['accuracy']*100:.1f}% "
                    f"(5-fold stratified cross-validation on {len(INC_GROUPS)}-group, "
                    f"8,077-plasmid training set)")

                # Per-Inc-group table
                if cv_data["per_class"]:
                    rows = []
                    for inc, m in sorted(cv_data["per_class"].items()):
                        rows.append({
                            "Inc Group": inc,
                            "Precision": f"{m['precision']*100:.1f}%",
                            "Recall": f"{m['recall']*100:.1f}%",
                            "F1 Score": f"{m['f1']*100:.1f}%",
                            "Training Samples": m["support"],
                        })
                    cv_df = pd.DataFrame(rows)
                    st.dataframe(cv_df, use_container_width=True, hide_index=True)

                # Confusion pair warnings for current batch
                records = st.session_state.get("records", [])
                confusion_count = sum(1 for r in records if r.get("inc_confusion_pair"))
                if confusion_count > 0:
                    st.warning(
                        f"**{confusion_count} classification(s) involve known confusion pairs** — "
                        f"the top-2 Inc candidates for these plasmids are known to be frequently "
                        f"confused by the KNN classifier. Consider these assignments with lower confidence.")
                    confused = [{"Plasmid": r["plasmid_id"],
                                 "Predicted": r.get("inc_best_match", "?"),
                                 "Note": r.get("inc_confusion_note", "")}
                                for r in records if r.get("inc_confusion_pair")]
                    st.dataframe(pd.DataFrame(confused), use_container_width=True, hide_index=True)

                # Summary of confusion pairs in training
                if cv_data["confusion_pairs"]:
                    st.markdown(f"**{len(cv_data['confusion_pairs'])} known confusion pairs** "
                                f"(>5% misclassification rate between pair):")
                    top_pairs = sorted(cv_data["confusion_pairs"],
                                       key=lambda x: x["errors"], reverse=True)[:10]
                    cp_rows = [{"Pair": f"{cp['inc_a']} ↔ {cp['inc_b']}",
                                "CV Errors": cp["errors"],
                                "% of A": f"{cp['pct_a']:.1f}%",
                                "% of B": f"{cp['pct_b']:.1f}%"}
                               for cp in top_pairs]
                    st.dataframe(pd.DataFrame(cp_rows), use_container_width=True, hide_index=True)

        # ── L3: Assembly Completeness ─────────────────────────────────────
        completeness_df = st.session_state.get("completeness_df")
        if completeness_df is not None and len(completeness_df) > 0:
            with st.expander("Assembly Completeness Assessment"):
                n_complete = (completeness_df["completeness_status"] == "COMPLETE").sum()
                n_near = (completeness_df["completeness_status"] == "NEAR-COMPLETE").sum()
                n_frag = (completeness_df["completeness_status"] == "FRAGMENTED").sum()
                n_poor = (completeness_df["completeness_status"] == "POOR").sum()
                c1, c2, c3, c4 = st.columns(4)
                c1.metric("Complete", n_complete)
                c2.metric("Near-complete", n_near)
                c3.metric("Fragmented", n_frag)
                c4.metric("Poor", n_poor)
                if n_frag + n_poor > 0:
                    st.warning(f"{n_frag + n_poor} plasmid(s) may have fragmented assemblies. "
                               "4-mer classification may be less reliable for these.")
                st.dataframe(completeness_df, use_container_width=True, hide_index=True)

        # ── L4: Database Coverage ─────────────────────────────────────────
        coverage_df = st.session_state.get("coverage_df")
        if coverage_df is not None and len(coverage_df) > 0:
            with st.expander("Database Coverage & Novelty Assessment"):
                n_green = (coverage_df["coverage_indicator"] == "GREEN").sum()
                n_yellow = (coverage_df["coverage_indicator"] == "YELLOW").sum()
                n_red = (coverage_df["coverage_indicator"] == "RED").sum()
                c1, c2, c3 = st.columns(3)
                c1.metric("Well-covered", n_green)
                c2.metric("Sparse coverage", n_yellow)
                c3.metric("Potentially novel", n_red)
                if n_red > 0:
                    st.error(f"{n_red} plasmid(s) may represent lineages not well-covered "
                             "by the reference database.")
                st.dataframe(coverage_df, use_container_width=True, hide_index=True)

        # ── L6: Novel Inc Group Discovery ─────────────────────────────────
        novel_df = st.session_state.get("novel_inc_df")
        if novel_df is not None and len(novel_df) > 0:
            with st.expander("Novel Inc/Rep Group Discovery"):
                putative = novel_df[novel_df["is_putative_novel_group"]]
                if len(putative) > 0:
                    n_groups = putative["novel_cluster"].nunique()
                    st.success(f"Discovered {n_groups} putative novel Inc/Rep group(s) "
                               f"from {len(putative)} plasmids with low-confidence classification.")
                else:
                    st.info("No putative novel groups found (requires >=3 unclassified plasmids "
                            "clustering tightly together).")
                st.dataframe(novel_df, use_container_width=True, hide_index=True)

        # ── L10: Gene Architecture / MGE Boundaries ───────────────────────
        mge_df = st.session_state.get("mge_df")
        if mge_df is not None and len(mge_df) > 0:
            with st.expander("Gene Architecture & MGE Boundaries"):
                # Summary
                n_is = (mge_df["mge_subtype"] == "IS_element").sum()
                n_int = (mge_df["mge_subtype"] == "integrase").sum()
                n_res = (mge_df["gene_type"] == "resistance").sum()
                n_backbone = (mge_df["gene_type"] == "backbone").sum()
                c1, c2, c3, c4 = st.columns(4)
                c1.metric("IS elements", n_is)
                c2.metric("Integrases", n_int)
                c3.metric("Resistance genes", n_res)
                c4.metric("Backbone genes", n_backbone)

                # Gene map for first plasmid
                records = st.session_state.get("records", [])
                source_files = mge_df["source_file"].unique()
                if len(source_files) > 0 and len(records) > 0:
                    selected_plasmid = st.selectbox(
                        "Select plasmid for gene map:",
                        source_files, key="mge_map_select"
                    )
                    plasmid_mge = mge_df[mge_df["source_file"] == selected_plasmid]
                    # Find plasmid length
                    p_len = 0
                    for r in records:
                        if selected_plasmid in r["plasmid_id"] or r["plasmid_id"] in selected_plasmid:
                            p_len = r["length"]
                            break
                    if p_len > 0 and len(plasmid_mge) > 0:
                        fig = draw_plasmid_gene_map(plasmid_mge, p_len, selected_plasmid)
                        st.pyplot(fig)
                        plt.close(fig)

                st.dataframe(
                    mge_df[["source_file", "gene_name", "start", "end", "strand",
                            "gene_type", "mge_subtype"]],
                    use_container_width=True, hide_index=True
                )


# ── TAB 3: Cladogram ────────────────────────────────────────────────────────

with tab_clado:
    if not st.session_state.analysis_done:
        st.info("Run analysis first to see cladograms.")
    elif st.session_state.Z is None:
        st.header("Cladogram Visualization")
        st.info(
            "Cladogram visualization requires 2 or more plasmids. "
            "In single-plasmid query mode, pLIN codes are assigned by "
            "nearest-neighbour lookup against the training database."
        )
        # Show query-mode summary instead
        qm = st.session_state.get("query_metadata")
        if qm:
            m = qm[0]
            st.markdown("### Query Mode Assignment Summary")
            col1, col2, col3 = st.columns(3)
            col1.metric("Nearest Neighbour", m["nn_plasmid"])
            col2.metric("Cosine Distance", f"{m['nn_distance']:.6f}")
            col3.metric("Neighbour Inc Type", m["nn_inc_type"])
            st.markdown(f"**Neighbour pLIN:** `{m['nn_plin']}`")
            st.markdown(f"**Assigned pLIN:** `{st.session_state.plin_codes[0]}`")

            # Show threshold breakdown
            thresholds = st.session_state.get("active_thresholds", PLIN_THRESHOLDS)
            rows = []
            nn_plin_parts = m["nn_plin"].split(".")
            assigned_parts = st.session_state.plin_codes[0].split(".")
            for idx, (level, thresh) in enumerate(thresholds.items()):
                inherited = m["nn_distance"] <= thresh
                rows.append({
                    "Level": f"L{idx+1} (Bin {level})",
                    "Threshold": f"d <= {thresh:.3f}",
                    "Distance": f"{m['nn_distance']:.6f}",
                    "Status": "Inherited" if inherited else "New branch",
                    "Neighbour ID": nn_plin_parts[idx] if idx < len(nn_plin_parts) else "?",
                    "Assigned ID": assigned_parts[idx] if idx < len(assigned_parts) else "?",
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
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

        # ── Within-group sequence alignment ─────────────────────────────────
        st.markdown("---")
        st.subheader("🔬 Within-Group Sequence Alignment")
        st.caption(
            "A shared pLIN code means two plasmids look similar by composition. "
            "This does not by itself say *which* regions are shared, or whether "
            "AMR content differs. Select two or more plasmids — typically ones "
            "sharing a pLIN code or prefix — to run pairwise BLAST and see "
            "coverage, identity, and shared-block coordinates directly."
        )

        blastn_path, makeblastdb_path = detect_blastn()
        if not blastn_path or not makeblastdb_path:
            st.warning(
                "BLAST+ not detected (blastn/makeblastdb). Within-group alignment "
                "is optional and requires BLAST+ to be installed and on PATH."
            )
        else:
            plin_lookup = dict(zip(labels, plin_codes))
            group_options = sorted(set(plin_codes))
            preselect_group = st.selectbox(
                "Filter by pLIN code (optional — narrows the plasmid list below)",
                ["(no filter)"] + group_options, index=0, key="align_group_filter",
            )
            candidate_ids = (
                [pid for pid in labels if plin_lookup.get(pid) == preselect_group]
                if preselect_group != "(no filter)" else labels
            )

            selected_ids = st.multiselect(
                "Plasmids to align (all pairwise combinations will be run)",
                options=candidate_ids,
                default=candidate_ids[:2] if len(candidate_ids) >= 2 else candidate_ids,
                key="align_selected_ids",
            )

            if len(selected_ids) >= 2:
                n_pairs = len(selected_ids) * (len(selected_ids) - 1) // 2
                run_align = st.button(f"Run alignment ({n_pairs} pair{'s' if n_pairs != 1 else ''})",
                                       key="run_within_group_alignment")
                if run_align:
                    with st.spinner(f"Running blastn on {n_pairs} plasmid pair(s)..."):
                        align_result = run_within_group_alignment(
                            st.session_state.records, selected_ids, blastn_binary=blastn_path)
                    st.session_state.within_group_alignment = align_result

                align_result = st.session_state.get("within_group_alignment")
                if align_result:
                    if align_result["error"]:
                        st.error(align_result["error"])
                    elif len(align_result["pairs"]) > 0:
                        st.markdown("**Pairwise summary**")
                        pairs_display = align_result["pairs"].rename(columns={
                            "plasmid_1": "Plasmid 1", "plasmid_2": "Plasmid 2",
                            "length_1_bp": "Length 1 (bp)", "length_2_bp": "Length 2 (bp)",
                            "coverage_pct": "Coverage (%)", "weighted_identity_pct": "Weighted Identity (%)",
                            "n_blocks": "Aligned Blocks",
                        }).drop(columns=["error"], errors="ignore")
                        st.dataframe(pairs_display, use_container_width=True, hide_index=True)
                        st.download_button(
                            "📥 Download pairwise summary (TSV)",
                            align_result["pairs"].to_csv(sep="\t", index=False),
                            "within_group_alignment_pairs.tsv", "text/tab-separated-values",
                        )

                        amr_df = st.session_state.get("amr_df")
                        amr_annotated = (
                            annotate_alignment_blocks_with_amr(align_result["blocks"], amr_df)
                            if amr_df is not None and len(amr_df) > 0 else pd.DataFrame()
                        )

                        st.markdown("**Pair detail**")
                        pair_labels = [f"{r['plasmid_1']} vs {r['plasmid_2']}"
                                       for _, r in align_result["pairs"].iterrows()]
                        chosen_pair = st.selectbox("Select a pair to view", pair_labels,
                                                    key="align_pair_detail_select")
                        p1, p2 = chosen_pair.split(" vs ")
                        pair_row = align_result["pairs"][
                            (align_result["pairs"]["plasmid_1"] == p1) &
                            (align_result["pairs"]["plasmid_2"] == p2)
                        ].iloc[0]
                        blocks_for_pair = align_result["blocks"][
                            (align_result["blocks"]["plasmid_1"] == p1) &
                            (align_result["blocks"]["plasmid_2"] == p2)
                        ]

                        if len(blocks_for_pair) == 0:
                            st.info("No aligned blocks above the length cutoff for this pair — "
                                    "these plasmids do not share substantial sequence despite the "
                                    "shared pLIN code prefix.")
                        else:
                            amr_for_pair = (
                                amr_annotated[(amr_annotated["plasmid_id"].isin([p1, p2]))]
                                if len(amr_annotated) > 0 else None
                            )
                            fig_align = draw_alignment_comparison(
                                blocks_for_pair, p1, p2,
                                int(pair_row["length_1_bp"]), int(pair_row["length_2_bp"]),
                                amr_annotated=amr_for_pair,
                            )
                            st.pyplot(fig_align, use_container_width=True)
                            st.download_button(
                                "📥 Download comparison (PNG)", fig_to_bytes(fig_align, "png"),
                                f"alignment_{p1}_vs_{p2}.png", "image/png",
                            )
                            plt.close(fig_align)

                            with st.expander(f"Aligned block coordinates ({len(blocks_for_pair)})"):
                                st.dataframe(
                                    blocks_for_pair[["q_start", "q_end", "s_start", "s_end",
                                                     "length_bp", "pct_identity"]].rename(columns={
                                        "q_start": f"{p1} start", "q_end": f"{p1} end",
                                        "s_start": f"{p2} start", "s_end": f"{p2} end",
                                        "length_bp": "Length (bp)", "pct_identity": "% Identity",
                                    }),
                                    use_container_width=True, hide_index=True,
                                )

                            if amr_for_pair is not None and len(amr_for_pair) > 0:
                                with st.expander("AMR genes: inside vs. outside shared blocks"):
                                    st.dataframe(
                                        amr_for_pair[["plasmid_id", "gene", "class", "start", "end",
                                                      "in_shared_block"]].rename(columns={
                                            "plasmid_id": "Plasmid", "gene": "Gene", "class": "Class",
                                            "start": "Start", "end": "End",
                                            "in_shared_block": "Inside shared block",
                                        }),
                                        use_container_width=True, hide_index=True,
                                    )
            else:
                st.info("Select at least 2 plasmids above to run alignment.")


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

                # Show prediction source
                if "mobility_source" in mob_df.columns:
                    sources = mob_df["mobility_source"].value_counts()
                    source_str = ", ".join(f"{src}: {cnt}" for src, cnt in sources.items())
                    st.caption(f"Prediction source: {source_str}")

                # Pie chart
                fig_mob = px.pie(
                    values=mob_counts.values, names=mob_counts.index,
                    title="Mobility Classification",
                    color=mob_counts.index,
                    color_discrete_map=mob_colors,
                )
                fig_mob.update_layout(height=350)
                st.plotly_chart(fig_mob, use_container_width=True)

                # Relaxase family distribution (when MOBsuite data available)
                if "relaxase_family" in mob_df.columns:
                    relaxase_data = mob_df[mob_df["relaxase_family"].astype(str).str.strip() != ""]
                    if len(relaxase_data) > 0:
                        rel_counts = relaxase_data["relaxase_family"].value_counts()
                        fig_rel = px.bar(
                            x=rel_counts.index, y=rel_counts.values,
                            title="Relaxase (MOB) Family Distribution",
                            labels={"x": "Relaxase Family", "y": "Count"},
                            color=rel_counts.index,
                        )
                        fig_rel.update_layout(height=300, showlegend=False)
                        st.plotly_chart(fig_rel, use_container_width=True)

                # MPF type distribution (when MOBsuite data available)
                if "mpf_type" in mob_df.columns:
                    mpf_data = mob_df[mob_df["mpf_type"].astype(str).str.strip() != ""]
                    if len(mpf_data) > 0:
                        mpf_counts = mpf_data["mpf_type"].value_counts()
                        fig_mpf = px.bar(
                            x=mpf_counts.index, y=mpf_counts.values,
                            title="Mating Pair Formation (MPF) Type Distribution",
                            labels={"x": "MPF Type", "y": "Count"},
                            color=mpf_counts.index,
                        )
                        fig_mpf.update_layout(height=300, showlegend=False)
                        st.plotly_chart(fig_mpf, use_container_width=True)

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
                st.info("No mobility data. Run analysis with AMRFinderPlus or MOBsuite enabled.")

        # ─── Outbreak Detection ───
        with epi_col2:
            st.subheader("Outbreak / Clone Detection")
            outbreak_clusters = st.session_state.get("outbreak_clusters", [])

            if outbreak_clusters:
                st.metric("Suspected Outbreak Clusters", len(outbreak_clusters))

                risk_icons = {"HIGH": "🔴", "MODERATE": "🟡", "INFO": "🔵"}
                for i, cluster in enumerate(outbreak_clusters):
                    risk_color = risk_icons.get(cluster["risk_level"], "🔵")
                    with st.expander(
                        f"{risk_color} Cluster {i+1}: pLIN {cluster['pLIN']} "
                        f"({cluster['n_plasmids']} plasmids, {cluster['n_amr_genes']} AMR genes)"
                    ):
                        st.markdown(f"**Risk level:** {cluster['risk_level']}")
                        st.markdown(f"**Shared pLIN code:** {cluster['pLIN']}")
                        if cluster["amr_genes"]:
                            consistency = "identical across all plasmids" if cluster["amr_profile_consistent"] \
                                else "genes pooled — profiles differ slightly between plasmids"
                            st.markdown(f"**Shared AMR genes ({consistency}):** {', '.join(cluster['amr_genes'])}")
                        else:
                            st.markdown("**Shared AMR genes:** none detected / AMRFinderPlus not run")
                        st.markdown(f"**Plasmids:**")
                        for p in cluster["plasmids"]:
                            st.markdown(f"  - `{p}`")

                if any(c["risk_level"] == "HIGH" for c in outbreak_clusters):
                    st.error(
                        "**HIGH-RISK outbreak cluster(s) detected.** "
                        "Plasmids sharing identical pLIN codes AND "
                        "the same AMR resistance profile may indicate clonal spread."
                    )
                elif any(c["risk_level"] == "INFO" for c in outbreak_clusters):
                    st.info(
                        "Plasmids sharing an identical pLIN code were found (likely the same plasmid), "
                        "but no AMR gene data is available to corroborate resistance profile. "
                        "Run AMRFinderPlus for stronger evidence."
                    )
            else:
                st.info(
                    "No outbreak clusters detected. Outbreak detection flags groups of "
                    "plasmids sharing the same pLIN L6 code (F-level) AND identical "
                    "AMR resistance profiles."
                )

        # ─── Risk Summary ───
        st.divider()
        st.subheader("Dissemination Risk Summary")

        integrated = st.session_state.integrated_df
        mob_df = st.session_state.get("mobility_results")

        if integrated is not None and mob_df is not None and len(mob_df) > 0:
            # Build merge columns list based on available data
            mob_merge_cols = ["plasmid_id", "mobility"]
            if "relaxase_family" in mob_df.columns:
                mob_merge_cols.append("relaxase_family")
            if "mpf_type" in mob_df.columns:
                mob_merge_cols.append("mpf_type")

            risk_df = integrated[["plasmid_id", "pLIN", "inc_type", "AMR_count"]].merge(
                mob_df[mob_merge_cols], on="plasmid_id", how="left"
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
                # Show columns dynamically based on available data
                display_cols = ["plasmid_id", "pLIN", "inc_type", "AMR_count", "mobility"]
                if "relaxase_family" in high_risk.columns:
                    display_cols.append("relaxase_family")
                if "mpf_type" in high_risk.columns:
                    display_cols.append("mpf_type")
                st.dataframe(high_risk[display_cols],
                             use_container_width=True, hide_index=True)

        # ─── Temporal Outbreak Clusters ───
        temporal_clusters = st.session_state.get("temporal_outbreak_clusters", [])
        if temporal_clusters:
            st.divider()
            st.subheader("Temporal Outbreak Clusters")
            st.caption("Plasmids sharing an identical pLIN code with collection dates within a 30-day window (matching AMR genes shown as supporting evidence)")
            st.metric("Temporal Outbreak Clusters", len(temporal_clusters))

            for i, tc in enumerate(temporal_clusters):
                risk_icon = {"CRITICAL": "🔴", "HIGH": "🟠", "MODERATE": "🟡"}.get(tc["risk_level"], "⚪")
                with st.expander(
                    f"{risk_icon} Cluster {i+1}: pLIN {tc['pLIN']} — "
                    f"{tc['n_plasmids']} plasmids, {tc['date_range_days']}d span, "
                    f"{tc['n_amr_genes']} AMR genes"
                ):
                    st.markdown(f"**Risk level:** {tc['risk_level']}")
                    st.markdown(f"**Date range:** {tc['earliest_date']} → {tc['latest_date']} ({tc['date_range_days']} days)")
                    if tc.get("locations"):
                        st.markdown(f"**Locations:** {', '.join(str(l) for l in tc['locations'])}")
                    st.markdown(f"**Shared AMR genes:** {', '.join(tc['amr_genes']) if tc['amr_genes'] else 'None'}")
                    st.markdown("**Plasmids:**")
                    for p in tc["plasmids"]:
                        st.markdown(f"  - `{p}`")

            if any(tc["risk_level"] == "CRITICAL" for tc in temporal_clusters):
                st.error(
                    "**CRITICAL temporal outbreak cluster(s) detected.** "
                    "Plasmids with identical L6 codes, same AMR profile, and collection dates "
                    "within 7 days — strongly suggestive of active clonal transmission."
                )

        # ─── Pathogen-Plasmid Integration (MLST + pLIN) ───
        mlst_df = st.session_state.get("mlst_df")
        transmission_pairs = st.session_state.get("transmission_pairs_df")
        transmission_summary = st.session_state.get("transmission_summary")

        if mlst_df is not None and len(mlst_df) > 0:
            st.divider()
            st.subheader("Pathogen-Plasmid Integration (MLST + pLIN)")
            st.caption(
                "Combined chromosomal typing (MLST) and plasmid typing (pLIN) "
                "to distinguish clonal vs horizontal plasmid spread"
            )

            with st.expander(f"MLST Typing Results ({len(mlst_df)} genomes)", expanded=True):
                mc1, mc2, mc3 = st.columns(3)
                mc1.metric("Genomes Typed", len(mlst_df))
                valid_st = mlst_df[mlst_df["ST"].astype(str).str.match(r"^\d+$")]
                mc2.metric("Unique STs", valid_st["ST"].nunique() if len(valid_st) > 0 else 0)
                mc3.metric("Schemes Detected", mlst_df[mlst_df["scheme"] != "-"]["scheme"].nunique())
                st.dataframe(mlst_df, use_container_width=True, hide_index=True)

            if transmission_pairs is not None and len(transmission_pairs) > 0 and transmission_summary:
                with st.expander("Transmission Mode Analysis", expanded=True):
                    tc1, tc2, tc3, tc4 = st.columns(4)
                    tc1.metric(
                        "Clonal Spread",
                        transmission_summary.get("clonal_spread", 0),
                        help="Same ST + same pLIN = vertical transmission",
                    )
                    tc2.metric(
                        "Horizontal Transfer",
                        transmission_summary.get("horizontal_transfer", 0),
                        help="Different STs + same pLIN = HGT",
                    )
                    tc3.metric(
                        "Same Strain, Diff Plasmid",
                        transmission_summary.get("same_strain_diff_plasmid", 0),
                    )
                    tc4.metric(
                        "Independent",
                        transmission_summary.get("independent", 0),
                    )

                    # Pie chart
                    mode_counts = transmission_pairs["transmission_mode"].value_counts()
                    if len(mode_counts) > 0:
                        mode_colors = {
                            "Clonal spread": "#E53935",
                            "Horizontal plasmid transfer": "#FB8C00",
                            "Same strain, different plasmids": "#1E88E5",
                            "Independent": "#43A047",
                        }
                        fig_trans = px.pie(
                            values=mode_counts.values,
                            names=mode_counts.index,
                            title="Transmission Mode Distribution",
                            color=mode_counts.index,
                            color_discrete_map=mode_colors,
                        )
                        fig_trans.update_layout(height=350)
                        st.plotly_chart(fig_trans, use_container_width=True)

                    st.dataframe(
                        transmission_pairs[[
                            "genome_A", "genome_B", "ST_A", "ST_B",
                            "plasmid_A", "plasmid_B", "pLIN_A", "pLIN_B",
                            "transmission_mode",
                        ]],
                        use_container_width=True,
                        hide_index=True,
                    )

                    hgt_count = transmission_summary.get("horizontal_transfer", 0)
                    if hgt_count > 0:
                        st.warning(
                            f"**{hgt_count} horizontal plasmid transfer event(s) detected.** "
                            "Different bacterial strains carry the same plasmid (same pLIN L6), "
                            "indicating active plasmid dissemination requiring enhanced infection "
                            "control measures beyond standard contact precautions."
                        )
                    clonal_count = transmission_summary.get("clonal_spread", 0)
                    if clonal_count > 0:
                        st.error(
                            f"**{clonal_count} clonal spread event(s) detected.** "
                            "Same bacterial strain (MLST ST) carrying the same plasmid (pLIN), "
                            "indicating direct person-to-person transmission."
                        )

        # ─── ANI Validation Results ───
        mash_df = st.session_state.get("mash_df")
        fastani_df = st.session_state.get("fastani_df")
        snp_df = st.session_state.get("snp_subtype_df")

        if (mash_df is not None and len(mash_df) > 0) or \
           (fastani_df is not None and len(fastani_df) > 0) or \
           (snp_df is not None and len(snp_df) > 0):
            st.divider()
            st.subheader("ANI Validation & SNP Sub-typing")

        if mash_df is not None and len(mash_df) > 0:
            with st.expander(f"Mash ANI Estimates ({len(mash_df)} pairs)"):
                mc1, mc2, mc3 = st.columns(3)
                mc1.metric("Mean ANI", f"{mash_df['ani_estimate'].mean():.1f}%")
                mc2.metric("Min ANI", f"{mash_df['ani_estimate'].min():.1f}%")
                mc3.metric("Max ANI", f"{mash_df['ani_estimate'].max():.1f}%")
                st.dataframe(
                    mash_df[["query", "reference", "mash_distance", "ani_estimate", "p_value"]].sort_values("ani_estimate", ascending=False),
                    use_container_width=True, hide_index=True,
                )

        # ANI Concordance Check (cross-validates cosine distance vs Mash ANI)
        concordance_df = st.session_state.get("concordance_df")
        if concordance_df is not None and len(concordance_df) > 0:
            n_concordant = (concordance_df["concordance"] == "concordant").sum()
            n_discordant = (concordance_df["concordance"] == "discordant").sum()
            n_total = n_concordant + n_discordant
            if n_total > 0:
                with st.expander(f"Cosine-ANI Concordance Check ({n_concordant}/{n_total} concordant)"):
                    if n_discordant > 0:
                        st.warning(
                            f"**{n_discordant} discordant assignment(s)** — cosine distance and Mash ANI "
                            f"disagree for these plasmids. The 4-mer composition proxy may be unreliable "
                            f"for these specific sequences. Verify with FastANI or minimap2 alignment.")
                        disc = concordance_df[concordance_df["concordance"] == "discordant"]
                        st.dataframe(disc[["plasmid_id", "nn_distance", "mash_ani", "note"]],
                                     use_container_width=True, hide_index=True)
                    else:
                        st.success(
                            f"All {n_total} pLIN assignments are concordant with Mash ANI estimates — "
                            f"4-mer cosine distances correlate with sequence-level similarity.")

        if fastani_df is not None and len(fastani_df) > 0:
            with st.expander(f"FastANI True ANI ({len(fastani_df)} pairs)"):
                fc1, fc2, fc3 = st.columns(3)
                fc1.metric("Mean ANI", f"{fastani_df['ani'].mean():.1f}%")
                fc2.metric("Min ANI", f"{fastani_df['ani'].min():.1f}%")
                fc3.metric("Max ANI", f"{fastani_df['ani'].max():.1f}%")
                st.dataframe(
                    fastani_df[["query", "reference", "ani", "orthologous_matches", "total_fragments"]].sort_values("ani", ascending=False),
                    use_container_width=True, hide_index=True,
                )

        if snp_df is not None and len(snp_df) > 0:
            with st.expander(f"SNP Sub-typing within L6 Clusters ({len(snp_df)} comparisons)"):
                valid_snps = snp_df[snp_df["snp_count"] >= 0]
                if len(valid_snps) > 0:
                    sc1, sc2, sc3 = st.columns(3)
                    sc1.metric("L6 Clusters Analyzed", valid_snps["l6_cluster"].nunique())
                    sc2.metric("Mean SNPs", f"{valid_snps['snp_count'].mean():.1f}")
                    sc3.metric("Max SNPs", int(valid_snps["snp_count"].max()))

                    # Flag identical plasmids (0 SNPs)
                    identical = valid_snps[valid_snps["snp_count"] == 0]
                    non_ref_identical = identical[identical["plasmid_id"] != identical["reference_id"]]
                    if len(non_ref_identical) > 0:
                        st.warning(
                            f"**{len(non_ref_identical)} plasmid pair(s) with 0 SNP differences** — "
                            "these are likely identical or near-identical sequences, strongly suggesting "
                            "recent clonal transmission or the same plasmid isolated multiple times."
                        )

                st.dataframe(
                    snp_df[["l6_cluster", "plasmid_id", "reference_id", "snp_count", "alignment_identity"]],
                    use_container_width=True, hide_index=True,
                )

        # ── L5: Recombination Signals ─────────────────────────────────────
        recom_df = st.session_state.get("recombination_df")
        if recom_df is not None and len(recom_df) > 0:
            with st.expander("Recombination Signals"):
                n_high = (recom_df["recombination_flag"] == "High").sum()
                n_med = (recom_df["recombination_flag"] == "Medium").sum()
                n_low = (recom_df["recombination_flag"] == "Low").sum()
                c1, c2, c3 = st.columns(3)
                c1.metric("High recombination", n_high)
                c2.metric("Medium", n_med)
                c3.metric("Low", n_low)
                if n_high > 0:
                    st.warning(f"{n_high} plasmid(s) show strong recombination signals "
                               "(low alignment coverage to nearest neighbor). "
                               "4-mer-based classification may be unreliable for these.")
                st.dataframe(recom_df, use_container_width=True, hide_index=True)

        # ── L7: Evolutionary Rate Estimation ──────────────────────────────
        evo_rate_df = st.session_state.get("evo_rate_df")
        if evo_rate_df is not None and len(evo_rate_df) > 0:
            with st.expander("Evolutionary Rate Estimation (Molecular Clock)"):
                st.markdown("SNP accumulation rates for L6 clusters with dated samples:")
                st.dataframe(evo_rate_df, use_container_width=True, hide_index=True)

        # ── L8: Cluster Stability Assessment ──────────────────────────────
        stability_df = st.session_state.get("stability_df")
        if stability_df is not None and len(stability_df) > 0:
            with st.expander("Cluster Robustness (Bootstrap Stability)"):
                n_stable = (stability_df["stability_status"] == "Stable").sum()
                n_moderate = (stability_df["stability_status"] == "Moderate").sum()
                n_unstable = (stability_df["stability_status"] == "Unstable").sum()
                c1, c2, c3 = st.columns(3)
                c1.metric("Stable (>80%)", n_stable)
                c2.metric("Moderate (50-80%)", n_moderate)
                c3.metric("Unstable (<50%)", n_unstable)
                if n_unstable > 0:
                    st.warning(f"{n_unstable} cluster(s) have low bootstrap support. "
                               "These assignments may change with additional data.")
                st.dataframe(stability_df, use_container_width=True, hide_index=True)

        linkage_df = st.session_state.get("linkage_comparison_df")
        if linkage_df is not None and len(linkage_df) > 0:
            with st.expander("Linkage Method Comparison"):
                st.markdown("Adjusted Rand Index (ARI) between clustering methods:")
                st.dataframe(linkage_df, use_container_width=True, hide_index=True)


# ── TAB 6: CRISPR Host Inference ────────────────────────────────────────────

with tab_crispr:
    if not st.session_state.analysis_done:
        st.info("Run analysis first to see CRISPR host inference results.")
    else:
        spacers_df = st.session_state.get("crispr_spacers_df")
        spacer_summary = st.session_state.get("crispr_spacer_summary_df")
        host_probs = st.session_state.get("crispr_host_probs_df")
        host_summary = st.session_state.get("crispr_host_summary_df")
        filtered_hits = st.session_state.get("crispr_filtered_hits_df")

        has_results = (
            host_probs is not None
            and isinstance(host_probs, pd.DataFrame)
            and len(host_probs) > 0
        )

        if not has_results:
            st.header("🧫 CRISPR Host Inference")
            st.info(
                "**No CRISPR host inference results available.**\n\n"
                "To enable CRISPR-based host prediction:\n"
                "1. Ensure **minced** and **BLAST** are installed\n"
                "2. Check **Run CRISPR host inference** in the sidebar\n"
                "3. Upload one or more **bacterial genome FASTA** files\n"
                "4. Re-run the analysis"
            )
        else:
            st.header("🧫 CRISPR Host Inference")
            st.caption(
                "Plasmid-host relationships inferred from CRISPR spacer matching "
                "(minced + BLASTN-short + softmax probability ranking)"
            )

            # ── Metrics row ──────────────────────────────────────────────
            n_hosts = spacer_summary["genome"].nunique() if spacer_summary is not None and len(spacer_summary) > 0 else 0
            n_arrays = spacer_summary["arrays"].sum() if spacer_summary is not None and len(spacer_summary) > 0 else 0
            n_spacers = len(spacers_df) if spacers_df is not None else 0
            n_predictions = len(host_probs)

            mc1, mc2, mc3, mc4 = st.columns(4)
            mc1.metric("Host Genomes Screened", n_hosts)
            mc2.metric("CRISPR Arrays Found", int(n_arrays))
            mc3.metric("Total Spacers", n_spacers)
            mc4.metric("Host–Plasmid Predictions", n_predictions)

            # ── High-confidence alert ────────────────────────────────────
            if "confidence" in host_probs.columns:
                high_conf = host_probs[host_probs["confidence"] == "High"]
                if len(high_conf) > 0:
                    st.success(
                        f"**{len(high_conf)} high-confidence host–plasmid prediction(s)** detected "
                        f"(probability >= 0.7). These represent strong CRISPR-based evidence of "
                        f"plasmid-host association."
                    )

            # ── Two-column layout ────────────────────────────────────────
            col_left, col_right = st.columns(2)

            with col_left:
                # Confidence distribution pie chart
                if "confidence" in host_probs.columns:
                    st.subheader("Prediction Confidence")
                    conf_counts = host_probs["confidence"].value_counts().reset_index()
                    conf_counts.columns = ["Confidence", "Count"]
                    color_map = {"High": "#2ecc71", "Medium": "#f39c12", "Low": "#e74c3c"}
                    fig_pie = px.pie(
                        conf_counts, names="Confidence", values="Count",
                        color="Confidence", color_discrete_map=color_map,
                        title="Confidence Distribution"
                    )
                    fig_pie.update_layout(height=350)
                    st.plotly_chart(fig_pie, use_container_width=True)

                # Top host predictions table
                st.subheader("Top Host Predictions")
                display_cols_hp = ["host_genome", "plasmid_id", "probability", "unique_spacers"]
                if "confidence" in host_probs.columns:
                    display_cols_hp.append("confidence")
                top_preds = host_probs.nlargest(20, "probability")
                st.dataframe(
                    top_preds[display_cols_hp],
                    use_container_width=True, hide_index=True
                )

            with col_right:
                # Spacers per host bar chart
                if spacer_summary is not None and len(spacer_summary) > 0:
                    st.subheader("Spacers per Host Genome")
                    fig_bar = px.bar(
                        spacer_summary.sort_values("spacers", ascending=False).head(20),
                        x="genome", y="spacers",
                        color="arrays", color_continuous_scale="Viridis",
                        labels={"genome": "Host Genome", "spacers": "Spacers", "arrays": "Arrays"},
                        title="CRISPR Spacers Extracted per Host"
                    )
                    fig_bar.update_layout(height=350, xaxis_tickangle=-45)
                    st.plotly_chart(fig_bar, use_container_width=True)

                # Host-plasmid probability heatmap
                st.subheader("Host–Plasmid Probability Heatmap")
                if len(host_probs) > 0:
                    # Pivot for heatmap (limit to top hosts/plasmids)
                    top_hosts = host_probs.groupby("host_genome")["probability"].max().nlargest(15).index
                    top_plasmids = host_probs.groupby("plasmid_id")["probability"].max().nlargest(20).index
                    hm_data = host_probs[
                        host_probs["host_genome"].isin(top_hosts)
                        & host_probs["plasmid_id"].isin(top_plasmids)
                    ]
                    if len(hm_data) > 0:
                        hm_pivot = hm_data.pivot_table(
                            index="host_genome", columns="plasmid_id",
                            values="probability", fill_value=0
                        )
                        fig_hm = px.imshow(
                            hm_pivot, aspect="auto",
                            color_continuous_scale="YlOrRd",
                            labels=dict(x="Plasmid", y="Host Genome", color="Probability"),
                            title="Host–Plasmid Association Probabilities"
                        )
                        fig_hm.update_layout(height=450)
                        st.plotly_chart(fig_hm, use_container_width=True)
                    else:
                        st.info("Not enough data for heatmap visualization.")

            # ── Expanders ────────────────────────────────────────────────
            with st.expander("Filtered BLAST Hits", expanded=False):
                if filtered_hits is not None and len(filtered_hits) > 0:
                    st.dataframe(filtered_hits, use_container_width=True, hide_index=True)
                else:
                    st.info("No filtered BLAST hits available.")

            with st.expander("All Extracted Spacers", expanded=False):
                if spacers_df is not None and len(spacers_df) > 0:
                    st.dataframe(spacers_df, use_container_width=True, hide_index=True)
                else:
                    st.info("No spacers extracted.")

            with st.expander("Full Probability Rankings", expanded=False):
                if host_probs is not None and len(host_probs) > 0:
                    st.dataframe(
                        host_probs.sort_values("probability", ascending=False),
                        use_container_width=True, hide_index=True
                    )
                else:
                    st.info("No probability rankings available.")


# ── TAB 7: DRAGNOME Buddy ───────────────────────────────────────────────────

with tab_buddy:
    st.header("🧬 DRAGNOME Buddy")
    st.caption("Your AI assistant for plasmid biology and AMR analysis")

    if not ollama_available:
        st.warning(
            "**Ollama not detected.** DRAGNOME Buddy requires Ollama running locally.\n\n"
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
            st.success("Analysis data loaded. DRAGNOME Buddy can answer questions about your results!", icon="✅")
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
            if prompt := st.chat_input("Ask DRAGNOME Buddy about your plasmids..."):
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


# ── TAB 8: Export ────────────────────────────────────────────────────────────

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

            # CRISPR host inference results
            crispr_host_summary = st.session_state.get("crispr_host_summary_df")
            crispr_spacers = st.session_state.get("crispr_spacers_df")
            crispr_host_probs = st.session_state.get("crispr_host_probs_df")
            if crispr_host_probs is not None and len(crispr_host_probs) > 0:
                crispr_probs_csv = crispr_host_probs.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 CRISPR Host Predictions (TSV)", crispr_probs_csv,
                                   "crispr_host_predictions.tsv", "text/tab-separated-values")
            if crispr_spacers is not None and len(crispr_spacers) > 0:
                crispr_spacers_csv = crispr_spacers.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 CRISPR Spacers (TSV)", crispr_spacers_csv,
                                   "crispr_spacers.tsv", "text/tab-separated-values")
            if crispr_host_summary is not None and len(crispr_host_summary) > 0:
                crispr_summ_csv = crispr_host_summary.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 CRISPR Host Summary (TSV)", crispr_summ_csv,
                                   "crispr_host_summary.tsv", "text/tab-separated-values")

            # MLST & Transmission Mode exports
            mlst_results = st.session_state.get("mlst_df")
            if mlst_results is not None and len(mlst_results) > 0:
                mlst_csv = mlst_results.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 MLST Results (TSV)", mlst_csv,
                                   "mlst_results.tsv", "text/tab-separated-values")
            trans_pairs = st.session_state.get("transmission_pairs_df")
            if trans_pairs is not None and len(trans_pairs) > 0:
                tp_csv = trans_pairs.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Transmission Mode Analysis (TSV)", tp_csv,
                                   "transmission_mode_analysis.tsv", "text/tab-separated-values")

            # Contig classification (plasmid vs chromosome)
            contig_cls = st.session_state.get("contig_classification")
            if contig_cls and len(contig_cls) > 0:
                contig_cls_df = pd.DataFrame(contig_cls)
                contig_csv = contig_cls_df.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Contig Classification (TSV)", contig_csv,
                                   "contig_classification.tsv", "text/tab-separated-values")

            # Excluded chromosomes
            excl_chr = st.session_state.get("excluded_chromosomes")
            if excl_chr is not None and len(excl_chr) > 0:
                excl_csv = excl_chr.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Excluded Chromosomes (TSV)", excl_csv,
                                   "excluded_chromosomes.tsv", "text/tab-separated-values")

            # Incomplete plasmids
            inc_plasm = st.session_state.get("incomplete_plasmids")
            if inc_plasm is not None and len(inc_plasm) > 0:
                inc_csv = inc_plasm.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Incomplete Plasmids (TSV)", inc_csv,
                                   "incomplete_plasmids.tsv", "text/tab-separated-values")

            # Assembly completeness
            compl_df = st.session_state.get("completeness_df")
            if compl_df is not None and len(compl_df) > 0:
                compl_csv = compl_df.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Assembly Completeness (TSV)", compl_csv,
                                   "assembly_completeness.tsv", "text/tab-separated-values")

            # Database coverage
            cov_df = st.session_state.get("coverage_df")
            if cov_df is not None and len(cov_df) > 0:
                cov_csv = cov_df.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Database Coverage (TSV)", cov_csv,
                                   "database_coverage.tsv", "text/tab-separated-values")

            # Novel Inc groups
            novel_df = st.session_state.get("novel_inc_df")
            if novel_df is not None and len(novel_df) > 0:
                novel_csv = novel_df.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Novel Inc Groups (TSV)", novel_csv,
                                   "novel_inc_groups.tsv", "text/tab-separated-values")

            # MGE boundaries
            mge_exp_df = st.session_state.get("mge_df")
            if mge_exp_df is not None and len(mge_exp_df) > 0:
                mge_csv = mge_exp_df.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 MGE Boundaries (TSV)", mge_csv,
                                   "mge_boundaries.tsv", "text/tab-separated-values")

            # Recombination signals
            recom_exp_df = st.session_state.get("recombination_df")
            if recom_exp_df is not None and len(recom_exp_df) > 0:
                recom_csv = recom_exp_df.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Recombination Signals (TSV)", recom_csv,
                                   "recombination_signals.tsv", "text/tab-separated-values")

            # Evolutionary rate
            evo_exp_df = st.session_state.get("evo_rate_df")
            if evo_exp_df is not None and len(evo_exp_df) > 0:
                evo_csv = evo_exp_df.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Evolutionary Rate (TSV)", evo_csv,
                                   "evolutionary_rate.tsv", "text/tab-separated-values")

            # Cluster stability
            stab_exp_df = st.session_state.get("stability_df")
            if stab_exp_df is not None and len(stab_exp_df) > 0:
                stab_csv = stab_exp_df.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Cluster Stability (TSV)", stab_csv,
                                   "cluster_stability.tsv", "text/tab-separated-values")

            # Linkage comparison
            link_exp_df = st.session_state.get("linkage_comparison_df")
            if link_exp_df is not None and len(link_exp_df) > 0:
                link_csv = link_exp_df.to_csv(sep="\t", index=False).encode()
                st.download_button("📥 Linkage Comparison (TSV)", link_csv,
                                   "linkage_comparison.tsv", "text/tab-separated-values")

        with col2:
            st.subheader("Figures")
            Z = st.session_state.Z
            labels = st.session_state.labels
            plin_codes = st.session_state.plin_codes
            strain_clusters = st.session_state.strain_clusters

            if Z is not None:
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
            else:
                st.info("Cladogram figures require 2+ plasmids.")

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
                    # CRISPR results
                    crispr_hp = st.session_state.get("crispr_host_probs_df")
                    crispr_sp = st.session_state.get("crispr_spacers_df")
                    crispr_hs = st.session_state.get("crispr_host_summary_df")
                    if crispr_hp is not None and len(crispr_hp) > 0:
                        zf.writestr("crispr_host_predictions.tsv",
                                    crispr_hp.to_csv(sep="\t", index=False))
                    if crispr_sp is not None and len(crispr_sp) > 0:
                        zf.writestr("crispr_spacers.tsv",
                                    crispr_sp.to_csv(sep="\t", index=False))
                    if crispr_hs is not None and len(crispr_hs) > 0:
                        zf.writestr("crispr_host_summary.tsv",
                                    crispr_hs.to_csv(sep="\t", index=False))
                    # MLST & Transmission Mode
                    mlst_res = st.session_state.get("mlst_df")
                    if mlst_res is not None and len(mlst_res) > 0:
                        zf.writestr("mlst_results.tsv",
                                    mlst_res.to_csv(sep="\t", index=False))
                    trans_p = st.session_state.get("transmission_pairs_df")
                    if trans_p is not None and len(trans_p) > 0:
                        zf.writestr("transmission_mode_analysis.tsv",
                                    trans_p.to_csv(sep="\t", index=False))
                    # Contig classification export
                    contig_cls = st.session_state.get("contig_classification")
                    if contig_cls and len(contig_cls) > 0:
                        zf.writestr("contig_classification.tsv",
                                    pd.DataFrame(contig_cls).to_csv(sep="\t", index=False))
                    excl_chr = st.session_state.get("excluded_chromosomes")
                    if excl_chr is not None and len(excl_chr) > 0:
                        zf.writestr("excluded_chromosomes.tsv",
                                    excl_chr.to_csv(sep="\t", index=False))
                    inc_plasm = st.session_state.get("incomplete_plasmids")
                    if inc_plasm is not None and len(inc_plasm) > 0:
                        zf.writestr("incomplete_plasmids.tsv",
                                    inc_plasm.to_csv(sep="\t", index=False))
                    # New analyses exports
                    for key, fname in [
                        ("completeness_df", "assembly_completeness.tsv"),
                        ("coverage_df", "database_coverage.tsv"),
                        ("novel_inc_df", "novel_inc_groups.tsv"),
                        ("mge_df", "mge_boundaries.tsv"),
                        ("recombination_df", "recombination_signals.tsv"),
                        ("evo_rate_df", "evolutionary_rate.tsv"),
                        ("stability_df", "cluster_stability.tsv"),
                        ("linkage_comparison_df", "linkage_comparison.tsv"),
                    ]:
                        zdf = st.session_state.get(key)
                        if zdf is not None and len(zdf) > 0:
                            zf.writestr(fname, zdf.to_csv(sep="\t", index=False))
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
