#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
FastANI validation of cosine-distance-based pLIN grouping, across ALL 28
Inc/Rep replicon groups (20 Gram-negative + 4 Gram-positive + 2 Acinetobacter
+ 2 Pseudomonas).

This reconstructs and generalises the one-off analysis that originally
produced source_data/Figure6_ANI_validation.tsv for the 20 Gram-negative
groups (npj AMR submission), so that a future re-run of the manuscript's
FastANI validation covers all 28 groups uniformly rather than only 20.

Method (matches the Gram-negative analysis, inferred from its output
structure: see manuscript Methods, "pLIN code assignment" / "Statistical
analysis"):
  1. Build a per-group sampling pool = training-set FASTAs
     (plasmid_sequences_for_training/<group>/fastas/) UNION high-confidence
     sequences from the classifier-predicted expanded reference database
     (output/reference_inc_classifications.tsv), restricted to a size band
     of [0.2x, 5x] the group's median sequence length (excludes fragments/
     outliers that FastANI cannot meaningfully align).
  2. Sample up to --sample-size (default 20) sequences per group
     (fixed --seed for reproducibility; groups with fewer members use all).
  3. Run fastANI (default params, one-to-one) on every directed pair
     (n*(n-1), excluding self).
  4. Compute the pipeline's own 4-mer cosine distance for the same pair
     (same method as assign_pLIN.py / assign_pLIN_reference.py).
  5. Pairs where FastANI returns no ANI value (below its internal
     detection floor) are recorded as such (ani_status=no_ANI_returned,
     fastani_ani left blank) rather than silently dropped or fabricated.

Usage:
  python validate_fastani_ani.py --groups all
  python validate_fastani_ani.py --groups repSA_large,repSA_small,repEF_conj,repEF_res,repAci1,repAci_large,repPae_large,repPae_small
  python validate_fastani_ani.py --groups IncFII,IncN --sample-size 20 --seed 42
  python validate_fastani_ani.py --fastani-binary /Users/basilxavier/miniforge3/bin/fastANI

Output:
  output/fastani_validation_<groups-label>.tsv   (per-pair results)
  output/fastani_validation_<groups-label>_summary.json (per-group + overall stats)
"""

import os
import sys
import glob
import json
import random
import shutil
import argparse
import subprocess
import tempfile
from itertools import product as iter_product

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import spearmanr

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TRAINING_DIR = os.path.join(BASE_DIR, "plasmid_sequences_for_training")
REFERENCE_DIR = os.path.join(BASE_DIR, "reference")
REF_CLASS_TSV = os.path.join(BASE_DIR, "output", "reference_inc_classifications.tsv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")

GRAM_NEGATIVE_GROUPS = [
    "ColE", "ColRNAI", "IncA", "IncAC2", "IncC", "IncF", "IncFIB", "IncFIBK",
    "IncFIC", "IncFII", "IncHI1", "IncHI2", "IncI", "IncI1", "IncI2", "IncN",
    "IncR", "IncX1", "IncX3", "IncX4",
]
GRAM_POSITIVE_NONFERMENTER_GROUPS = [
    "repSA_large", "repSA_small", "repEF_conj", "repEF_res",
    "repAci1", "repAci_large", "repPae_large", "repPae_small",
]
ALL_28_GROUPS = GRAM_NEGATIVE_GROUPS + GRAM_POSITIVE_NONFERMENTER_GROUPS

K = 4
BASES = "ACGT"
ALL_KMERS = ["".join(p) for p in iter_product(BASES, repeat=K)]
KMER_INDEX = {kmer: i for i, kmer in enumerate(ALL_KMERS)}
N_KMERS = len(ALL_KMERS)

DEFAULT_SAMPLE_SIZE = 20
DEFAULT_SEED = 42
SIZE_BAND_LOW = 0.2
SIZE_BAND_HIGH = 5.0
MIN_INC_CONFIDENCE = 0.9


# ── 4-mer cosine distance (matches assign_pLIN_reference.py:kmer_vector) ──────

def kmer_vector(sequence):
    """Compute normalised 4-mer frequency vector (sliding window)."""
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


def cosine_distance(v1, v2):
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denom == 0:
        return 1.0
    return 1.0 - float(np.dot(v1, v2) / denom)


# ── Sampling pool construction ─────────────────────────────────────────────────

def seq_length(fasta_path):
    try:
        rec = next(SeqIO.parse(fasta_path, "fasta"))
        return len(rec.seq)
    except StopIteration:
        return 0


def build_pool(group):
    """Sampling pool for one Inc/Rep group: training FASTAs UNION
    high-confidence expanded-reference sequences, size-band filtered."""
    raw_pool = {}

    train_dir = os.path.join(TRAINING_DIR, group, "fastas")
    for fp in sorted(glob.glob(os.path.join(train_dir, "*.fasta"))):
        acc = os.path.basename(fp)[:-len(".fasta")]
        raw_pool[acc] = fp

    # Optional staging dir of newly downloaded/verified expansion sequences
    # (see plasmid_sequences_for_training/<group>/fastas_expansion_new/),
    # unioned in the same way as the curated training-fastas loop above.
    expansion_dir = os.path.join(TRAINING_DIR, group, "fastas_expansion_new")
    if os.path.isdir(expansion_dir):
        for fp in sorted(glob.glob(os.path.join(expansion_dir, "*.fasta"))):
            acc = os.path.basename(fp)[:-len(".fasta")]
            if acc not in raw_pool:
                raw_pool[acc] = fp

    if os.path.exists(REF_CLASS_TSV):
        ref_df = pd.read_csv(REF_CLASS_TSV, sep="\t")
        sub = ref_df[
            (ref_df["inc_type"] == group)
            & (ref_df["inc_confidence"] >= MIN_INC_CONFIDENCE)
            & (ref_df["is_multiple"] == False)  # noqa: E712
        ]
        for acc in sub["plasmid_id"]:
            if acc in raw_pool:
                continue
            fp = os.path.join(REFERENCE_DIR, f"{acc}.fasta")
            if os.path.exists(fp):
                raw_pool[acc] = fp

    lengths = {acc: seq_length(fp) for acc, fp in raw_pool.items()}
    lens_sorted = sorted(lengths.values())
    if not lens_sorted:
        return {}
    median_len = lens_sorted[len(lens_sorted) // 2]
    lo, hi = SIZE_BAND_LOW * median_len, SIZE_BAND_HIGH * median_len
    return {acc: fp for acc, fp in raw_pool.items() if lo <= lengths[acc] <= hi}


def sample_group(group, pool, sample_size, seed):
    rng = random.Random(seed + hash(group) % 100000)
    ids = sorted(pool.keys())
    if len(ids) <= sample_size:
        return ids
    return sorted(rng.sample(ids, sample_size))


# ── FastANI ────────────────────────────────────────────────────────────────────

def detect_fastani_binary():
    """Same detection logic as plin_app.py:detect_fastani()."""
    which = shutil.which("fastANI")
    if which:
        return which
    for base, env in [
        (os.path.expanduser("~/miniforge3"), None),
        (os.path.expanduser("~/miniconda3/envs"), "pLIN_tools"),
    ]:
        candidate = os.path.join(base, env, "bin", "fastANI") if env else os.path.join(base, "bin", "fastANI")
        if os.path.exists(candidate):
            return candidate
    return None


def run_fastani_pair(fastani_binary, query_fp, ref_fp, tmpdir):
    out_path = os.path.join(tmpdir, "fastani_out.tsv")
    if os.path.exists(out_path):
        os.remove(out_path)
    cmd = [fastani_binary, "-q", query_fp, "-r", ref_fp, "-o", out_path]
    subprocess.run(cmd, capture_output=True, text=True)
    if not os.path.exists(out_path) or os.path.getsize(out_path) == 0:
        return None
    with open(out_path) as f:
        line = f.readline().strip()
    if not line:
        return None
    parts = line.split("\t")
    if len(parts) < 5:
        return None
    return {
        "ani": float(parts[2]),
        "ortho_matches": int(parts[3]),
        "total_frags": int(parts[4]),
    }


# ── Main pipeline ──────────────────────────────────────────────────────────────

def validate_groups(groups, sample_size, seed, fastani_binary, verbose=True):
    """Run the full pairwise cosine-distance vs FastANI-ANI validation for
    the given list of Inc/Rep groups. Returns (results_df, summary_dict)."""
    seq_cache = {}
    vec_cache = {}
    all_rows = []
    pool_info = {}

    with tempfile.TemporaryDirectory(prefix="fastani_validation_") as tmpdir:
        for group in groups:
            if verbose:
                print(f"=== {group} ===", flush=True)
            pool = build_pool(group)
            n_train = len(glob.glob(os.path.join(TRAINING_DIR, group, "fastas", "*.fasta")))
            pool_info[group] = {"pool_size_total": len(pool), "n_training_fastas": n_train}

            sample_ids = sample_group(group, pool, sample_size, seed)
            if verbose:
                print(f"  pool={len(pool)} sampled={len(sample_ids)}", flush=True)
            if len(sample_ids) < 2:
                pool_info[group].update(pairs_attempted=0, pairs_no_ani=0, pairs_with_ani=0)
                continue

            for acc in sample_ids:
                if acc in seq_cache:
                    continue
                fp = pool[acc]
                rec = next(SeqIO.parse(fp, "fasta"))
                seq = str(rec.seq)
                seq_cache[acc] = fp
                vec_cache[acc] = kmer_vector(seq)

            pairs = [(q, r) for q in sample_ids for r in sample_ids if q != r]
            n_no_ani = 0
            for q, r in pairs:
                fani = run_fastani_pair(fastani_binary, seq_cache[q], seq_cache[r], tmpdir)
                cdist = cosine_distance(vec_cache[q], vec_cache[r])
                row = {
                    "inc_type": group,
                    "query": q,
                    "reference": r,
                    "cosine_distance": cdist,
                }
                if fani is None:
                    n_no_ani += 1
                    row.update(fastani_ani=np.nan, ortho_matches=np.nan,
                               total_frags=np.nan, ani_status="no_ANI_returned")
                else:
                    row.update(fastani_ani=fani["ani"], ortho_matches=fani["ortho_matches"],
                               total_frags=fani["total_frags"], ani_status="ok")
                all_rows.append(row)

            pool_info[group].update(
                pairs_attempted=len(pairs),
                pairs_no_ani=n_no_ani,
                pairs_with_ani=len(pairs) - n_no_ani,
            )
            if verbose:
                print(f"  done: {len(pairs)-n_no_ani}/{len(pairs)} pairs with ANI value", flush=True)

    df = pd.DataFrame(all_rows)
    summary = _summarize(df, pool_info)
    return df, summary


def _summarize(df, pool_info):
    ok = df[df["ani_status"] == "ok"].copy() if len(df) else df
    per_group = []
    for g in sorted(df["inc_type"].unique()) if len(df) else []:
        sub = ok[ok["inc_type"] == g]
        n_attempt = len(df[df["inc_type"] == g])
        n_ok = len(sub)
        if n_ok >= 4:
            rho, p = spearmanr(sub["cosine_distance"], sub["fastani_ani"])
        else:
            rho, p = (None, None)
        per_group.append({
            "inc_type": g,
            **pool_info.get(g, {}),
            "yield_pct": round(100 * n_ok / n_attempt, 1) if n_attempt else None,
            "spearman_rho": round(rho, 4) if rho is not None else None,
            "spearman_p": p,
            "median_ani": round(sub["fastani_ani"].median(), 4) if n_ok else None,
            "min_ani": round(sub["fastani_ani"].min(), 4) if n_ok else None,
            "max_ani": round(sub["fastani_ani"].max(), 4) if n_ok else None,
        })

    overall = {}
    if len(ok):
        rho_all, p_all = spearmanr(ok["cosine_distance"], ok["fastani_ani"])
        l6 = ok[ok["cosine_distance"] <= 0.001]
        overall = {
            "pairs_attempted": int(len(df)),
            "pairs_with_ani_value": int(len(ok)),
            "yield_pct": round(100 * len(ok) / len(df), 1) if len(df) else None,
            "spearman_rho": round(rho_all, 4),
            "spearman_p": p_all,
            "median_ani_all_pairs": round(ok["fastani_ani"].median(), 4),
            "n_pairs_l6_equivalent": int(len(l6)),
            "median_ani_at_l6_equivalent": round(l6["fastani_ani"].median(), 4) if len(l6) else None,
        }
    return {"per_group": per_group, "overall": overall}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--groups", default="gram_positive_nonfermenter",
        help=(
            "Comma-separated Inc/Rep group names, or one of: "
            "'all' (all 28 groups), 'gram_negative' (20 groups), "
            "'gram_positive_nonfermenter' (8 groups, default)."
        ),
    )
    parser.add_argument("--sample-size", type=int, default=DEFAULT_SAMPLE_SIZE)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--fastani-binary", default=None, help="Path to fastANI binary (auto-detected if omitted).")
    parser.add_argument("--out-prefix", default=None, help="Output file prefix (default derived from --groups).")
    args = parser.parse_args()

    if args.groups == "all":
        groups = ALL_28_GROUPS
        label = "all28"
    elif args.groups == "gram_negative":
        groups = GRAM_NEGATIVE_GROUPS
        label = "gram_negative"
    elif args.groups == "gram_positive_nonfermenter":
        groups = GRAM_POSITIVE_NONFERMENTER_GROUPS
        label = "gram_positive_nonfermenter"
    else:
        groups = [g.strip() for g in args.groups.split(",") if g.strip()]
        label = "custom"

    fastani_binary = args.fastani_binary or detect_fastani_binary()
    if not fastani_binary:
        print("ERROR: fastANI binary not found. Install via `conda install -c bioconda fastani` "
              "or pass --fastani-binary /path/to/fastANI.", file=sys.stderr)
        sys.exit(1)
    print(f"Using fastANI binary: {fastani_binary}", flush=True)

    df, summary = validate_groups(groups, args.sample_size, args.seed, fastani_binary)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    prefix = args.out_prefix or f"fastani_validation_{label}"
    out_tsv = os.path.join(OUTPUT_DIR, f"{prefix}.tsv")
    out_json = os.path.join(OUTPUT_DIR, f"{prefix}_summary.json")
    df.to_csv(out_tsv, sep="\t", index=False)
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nWrote {len(df)} pair rows to {out_tsv}")
    print(f"Wrote summary to {out_json}")
    print(json.dumps(summary.get("overall", {}), indent=2, default=str))


if __name__ == "__main__":
    main()
