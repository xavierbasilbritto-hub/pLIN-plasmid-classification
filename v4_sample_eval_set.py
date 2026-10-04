#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Draw the pLIN v4 evaluation plasmids (output/backbone_v4/PREREGISTRATION.md, section 4).

2,000 plasmids uniformly at random (seed 2026) from the release database
(db-2026.10.02, every plasmid in output/pLIN_reference_assignments.tsv),
excluding the 2,747 plasmids of the earlier comparator benchmark. The same
random generator then shuffles the sample; the first 1,000 are the calibration
half, the rest the test half.

Usage:
  python v4_sample_eval_set.py
Output: output/backbone_v4/eval_plasmids.tsv (plasmid_id, fasta, half),
        output/backbone_v4/eval_{calibration,test,all}.txt (FASTA path lists)
"""

import os
import random

import pandas as pd

from validate_alignment_backbone import fasta_paths

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "output", "backbone_v4")
SEED, N = 2026, 2000


def main():
    release = sorted(set(pd.read_csv(os.path.join(BASE_DIR, "output", "pLIN_reference_assignments.tsv"),
                                     sep="\t", usecols=["plasmid_id"]).plasmid_id))
    training = fasta_paths()
    path = {i: training.get(i, os.path.join(BASE_DIR, "reference", f"{i}.fasta")) for i in release}
    missing = [i for i, p in path.items() if not os.path.exists(p)]
    assert not missing, f"{len(missing)} release plasmids have no FASTA, e.g. {missing[:3]}"

    earlier = {os.path.basename(l.strip())[:-len(".fasta")]
               for l in open(os.path.join(BASE_DIR, "output", "comparator_benchmark", "pling_inputs.txt")) if l.strip()}
    eligible = [i for i in release if i not in earlier]

    rng = random.Random(SEED)
    sample = rng.sample(eligible, N)
    rng.shuffle(sample)
    df = pd.DataFrame({"plasmid_id": sample, "fasta": [path[i] for i in sample],
                       "half": ["calibration"] * (N // 2) + ["test"] * (N - N // 2)})
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(os.path.join(OUT, "eval_plasmids.tsv"), sep="\t", index=False)
    for name, sub in (("calibration", df[df.half == "calibration"]), ("test", df[df.half == "test"]), ("all", df)):
        open(os.path.join(OUT, f"eval_{name}.txt"), "w").write("\n".join(sorted(sub.fasta)) + "\n")
    print(f"release {len(release)}, excluded {len(release) - len(eligible)}, eligible {len(eligible)}, "
          f"sampled {len(df)} ({(df.half == 'calibration').sum()} calibration / {(df.half == 'test').sum()} test)")


if __name__ == "__main__":
    main()
