#!/usr/bin/env python3
# Copyright (C) 2025-2026 Basil Britto Xavier. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Worked example: how the k-mer measures behave when two plasmids differ greatly in size.

Plasmid length spans four orders of magnitude, so a reviewer will reasonably ask how a measure
that is not size-invariant handles that. This takes a real release plasmid, extracts a segment of
it, and reports the asymmetric and symmetric similarities between the segment and the whole.

Usage:
  python v41_size_example.py
Output: output/backbone_v41/size_example.json
"""

import json
import os

from Bio import SeqIO

from plin_kmers import adaptive_sketch, containment, min_containment

BASE = os.path.dirname(os.path.abspath(__file__))
V41 = os.path.join(BASE, "output", "backbone_v41")
# a large release plasmid; the segment is a contiguous piece of it
SOURCE = "NZ_CP026661.1"
SEGMENT_BP = 10000


def main():
    path = os.path.join(BASE, "reference", f"{SOURCE}.fasta")
    seq = str(next(SeqIO.parse(path, "fasta")).seq)
    segment = seq[len(seq) // 4: len(seq) // 4 + SEGMENT_BP]
    a, sa = adaptive_sketch(segment)
    b, sb = adaptive_sketch(seq)
    out = {"source": SOURCE, "source_bp": len(seq), "segment_bp": len(segment),
           "containment_segment_in_whole": round(float(containment(a, b, sa, sb)), 3),
           "symmetric_similarity": round(float(min_containment(a, b, sa, sb)), 3),
           "segment_hashes": int(a.size), "whole_hashes": int(b.size)}
    json.dump(out, open(os.path.join(V41, "size_example.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
