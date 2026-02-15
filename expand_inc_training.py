#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Expand Inc-group training data by extracting labeled sequences from reference FASTA.

Parses FASTA headers from sequences.fasta to identify plasmid sequences with
Inc type labels, normalizes types, and copies qualifying sequences into the
training folder structure expected by build_inc_centroids.py.

Usage:  python expand_inc_training.py
"""

import os
import re
import sys
import shutil
from collections import defaultdict

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
REFERENCE_DIR = os.path.join(BASE_DIR, "reference")
TRAINING_DIR = os.path.join(BASE_DIR, "plasmid_sequences_for_training")
SEQUENCES_FASTA = os.path.join(BASE_DIR, "sequences.fasta")

# Minimum samples required per Inc group to include in training
MIN_SAMPLES = 10

# ── Canonical Inc type set ───────────────────────────────────────────────────
# Only these exact patterns are valid Inc types. This prevents false positives
# from species names like "incola", "incerta", "incanae", etc.

VALID_INC_TYPES = {
    # IncF family
    "IncF", "IncFI", "IncFIA", "IncFIB", "IncFIBK", "IncFIBpQil",
    "IncFIC", "IncFII", "IncFIIK",
    # IncH family
    "IncH", "IncHI1", "IncHI1B", "IncHI2", "IncHI2A",
    # IncI family
    "IncI", "IncI1", "IncI2",
    # IncX family
    "IncX", "IncX1", "IncX3", "IncX4",
    # IncA/C family
    "IncA", "IncAC2", "IncC",
    # Other well-characterized groups
    "IncB", "IncBOKZ", "IncK", "IncL", "IncLM",
    "IncN", "IncN3", "IncP", "IncQ", "IncQ1",
    "IncR", "IncU", "IncW", "IncY",
    # Col types
    "ColRNAI", "ColE", "ColE1", "Col440I", "ColKP3", "ColpVC",
}

# Normalization: map sub-variants to canonical training groups
NORMALIZE_MAP = {
    "IncFIBpQil": "IncFIB",
    "IncHI2A": "IncHI2",
    "IncHI1B": "IncHI1",
    "IncN3": "IncN",
    "IncQ1": "IncQ",
    "IncFI": "IncF",
    "IncH": "IncHI1",
    "IncLM": "IncL",
    "ColE1": "ColE",
    "Col440I": "ColE",
    "ColKP3": "ColRNAI",
    "ColpVC": "ColRNAI",
    "IncH12": "IncHI2",
}

# Regex to extract Inc type from FASTA header.
# Matches Inc types in plasmid name context (after "plasmid" keyword)
# Pattern: plasmid p<name>_IncXXX  or  plasmid p<name>-IncXXX  or  plasmid IncXXX
INC_PATTERN = re.compile(
    r'plasmid\s+\S*[_\-]?(Inc[A-Z][A-Za-z0-9]*|Col[A-Z][A-Za-z0-9]*)',
    re.IGNORECASE
)

# Standalone pattern for cases like "plasmid pIncI1_KP045"
INC_STANDALONE = re.compile(
    r'p(Inc[A-Z][A-Za-z0-9]*)[_\-]',
    re.IGNORECASE
)


def extract_inc_type(header):
    """Extract and validate Inc type from a FASTA header line.

    Returns normalized Inc type string or None if not found/invalid.
    """
    # Try main pattern first
    matches = INC_PATTERN.findall(header)
    if not matches:
        matches = INC_STANDALONE.findall(header)

    if not matches:
        return None

    # Take the first match
    raw_type = matches[0]

    # Normalize case: ensure "Inc" prefix is properly capitalized
    # Handle mixed case like "incFII", "incN", "inckii"
    if raw_type.lower().startswith("inc"):
        raw_type = "Inc" + raw_type[3:]
    elif raw_type.lower().startswith("col"):
        raw_type = "Col" + raw_type[3:]

    # Check against valid types (case-sensitive after normalization)
    if raw_type not in VALID_INC_TYPES:
        # Try common case fixes
        for valid in VALID_INC_TYPES:
            if raw_type.lower() == valid.lower():
                raw_type = valid
                break
        else:
            return None

    # Apply normalization map
    normalized = NORMALIZE_MAP.get(raw_type, raw_type)
    return normalized


def get_existing_accessions(inc_type):
    """Get set of accession IDs already in a training folder."""
    fastas_dir = os.path.join(TRAINING_DIR, inc_type, "fastas")
    if not os.path.isdir(fastas_dir):
        return set()
    accessions = set()
    for fname in os.listdir(fastas_dir):
        if fname.endswith(".fasta"):
            # Remove .fasta extension and any "RefSeq_" prefix
            acc = fname.replace(".fasta", "").replace("RefSeq_", "")
            accessions.add(acc)
    return accessions


def main():
    if not os.path.isdir(REFERENCE_DIR):
        print(f"ERROR: Reference directory not found: {REFERENCE_DIR}")
        print("Run the FASTA splitting step first.")
        sys.exit(1)

    # Parse headers from reference files (read only first line of each file)
    print("Scanning reference FASTA headers...")
    inc_assignments = defaultdict(list)  # {inc_type: [(accession, fasta_path), ...]}
    total_scanned = 0
    total_labeled = 0

    # Read headers from the original sequences.fasta for speed
    # (reading first line of 72k files would be slower)
    if os.path.exists(SEQUENCES_FASTA):
        print(f"Reading headers from {SEQUENCES_FASTA}...")
        with open(SEQUENCES_FASTA, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                total_scanned += 1
                if total_scanned % 10000 == 0:
                    print(f"  Scanned {total_scanned} headers...")

                header = line.strip()
                accession = header.split()[0][1:]  # Remove '>' prefix

                inc_type = extract_inc_type(header)
                if inc_type:
                    fasta_path = os.path.join(REFERENCE_DIR, f"{accession}.fasta")
                    if os.path.isfile(fasta_path):
                        inc_assignments[inc_type].append((accession, fasta_path))
                        total_labeled += 1
    else:
        print(f"sequences.fasta not found, scanning reference directory...")
        ref_files = sorted(os.listdir(REFERENCE_DIR))
        for fname in ref_files:
            if not fname.endswith(".fasta"):
                continue
            total_scanned += 1
            if total_scanned % 10000 == 0:
                print(f"  Scanned {total_scanned} files...")

            fpath = os.path.join(REFERENCE_DIR, fname)
            with open(fpath, "r") as f:
                header = f.readline().strip()

            accession = fname.replace(".fasta", "")
            inc_type = extract_inc_type(header)
            if inc_type:
                inc_assignments[inc_type].append((accession, fpath))
                total_labeled += 1

    print(f"\nScanned {total_scanned} sequences, found {total_labeled} with Inc type labels")
    print(f"Unique Inc types found: {len(inc_assignments)}")

    # Report all types and counts
    print("\n── Inc Type Counts (from reference) ──")
    for inc_type in sorted(inc_assignments.keys(), key=lambda x: -len(inc_assignments[x])):
        count = len(inc_assignments[inc_type])
        status = "OK" if count >= MIN_SAMPLES else f"SKIP (<{MIN_SAMPLES})"
        print(f"  {inc_type:15s} {count:5d}  {status}")

    # Copy qualifying sequences to training folders
    print(f"\n── Populating training folders (minimum {MIN_SAMPLES} samples) ──")
    groups_added = 0
    sequences_added = 0

    for inc_type in sorted(inc_assignments.keys()):
        entries = inc_assignments[inc_type]
        if len(entries) < MIN_SAMPLES:
            continue

        fastas_dir = os.path.join(TRAINING_DIR, inc_type, "fastas")
        os.makedirs(fastas_dir, exist_ok=True)

        # Get existing accessions to avoid duplicates
        existing = get_existing_accessions(inc_type)
        added = 0

        for accession, fasta_path in entries:
            # Check for duplicates (by accession, with and without version)
            acc_base = accession.split(".")[0]
            if accession in existing or acc_base in existing:
                continue

            dest = os.path.join(fastas_dir, f"{accession}.fasta")
            if not os.path.exists(dest):
                shutil.copy2(fasta_path, dest)
                added += 1

        existing_count = len(existing)
        total_count = existing_count + added
        print(f"  {inc_type:15s} existing: {existing_count:5d}  added: {added:4d}  total: {total_count:5d}")
        sequences_added += added
        groups_added += 1

    print(f"\n── Summary ──")
    print(f"  Inc groups with training data: {groups_added}")
    print(f"  New sequences added: {sequences_added}")
    print(f"  Training directory: {TRAINING_DIR}")
    print(f"\nNext step: run  python build_inc_centroids.py  to retrain the classifier.")


if __name__ == "__main__":
    main()
