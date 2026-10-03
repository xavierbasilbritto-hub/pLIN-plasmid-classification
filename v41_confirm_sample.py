#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Draw the pLIN v4.1 confirmatory data (PREREGISTRATION_v4.1.md, section 5).

Population: the 127,517 unique accessions of db-2026.10.03 (RefSeq_ duplicates
dropped, plain accession kept). Excluded (by accession): every development
plasmid: the 2,747 comparator-benchmark plasmids, the 2,000 v4 evaluation
plasmids, the 50 mutation-test plasmids (seed 3), the 200 re-query-check
plasmids (seed 11) and the outbreak plasmids (74 accessions).

Draws (each from what remains after the previous draw):
  test     2,000 plasmids, seed 4141 (pre-registered)
  mutation 100 plasmids,   seed 4242 (pre-registered; mutations use seed 4343)
  requery  500 plasmids,   seed 4444 (count pre-registered; seed fixed here,
           before any result)

Output: output/backbone_v41/confirm/{test,mutation,requery}_plasmids.tsv, test_all.txt
"""

import os
import random

import pandas as pd

import plin_v4 as v
from v41_sketch_all import unique_ids

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "output", "backbone_v41", "confirm")


def acc(i):
    return i.replace("RefSeq_", "", 1)


def main():
    os.makedirs(OUT, exist_ok=True)
    paths = v.release_plasmids()
    ids = unique_ids()
    dev = set()
    dev |= {os.path.basename(l.strip())[:-len(".fasta")]
            for l in open(os.path.join(BASE_DIR, "output", "comparator_benchmark", "pling_inputs.txt")) if l.strip()}
    dev |= set(pd.read_csv(os.path.join(BASE_DIR, "output", "backbone_v4", "eval_plasmids.tsv"), sep="\t").plasmid_id)
    allrel = sorted(paths)
    dev |= set(random.Random(3).sample([i for i in allrel if not i.startswith("RefSeq_")], 50))
    dev |= set(random.Random(11).sample(allrel, 200))
    ob = pd.read_csv(os.path.join(BASE_DIR, "output", "outbreak_validation_founder_results.tsv"), sep="\t")
    dev |= set(ob.accession) | set(ob.database_id.dropna())
    dev_acc = {acc(x) for x in dev}
    pool = [i for i in ids if acc(i) not in dev_acc]

    test = sorted(random.Random(4141).sample(pool, 2000))
    rest = [i for i in pool if i not in set(test)]
    mut = sorted(random.Random(4242).sample(rest, 100))
    rest = [i for i in rest if i not in set(mut)]
    req = sorted(random.Random(4444).sample(rest, 500))
    for name, sel in (("test", test), ("mutation", mut), ("requery", req)):
        pd.DataFrame({"plasmid_id": sel, "fasta": [paths[i] for i in sel]}) \
            .to_csv(os.path.join(OUT, f"{name}_plasmids.tsv"), sep="\t", index=False)
    open(os.path.join(OUT, "test_all.txt"), "w").write("\n".join(paths[i] for i in test) + "\n")
    print(f"population {len(ids):,}; development accessions excluded {len(ids) - len(pool):,}; "
          f"drawn: test {len(test)}, mutation {len(mut)}, requery {len(req)}")
    assert not (set(map(acc, test)) & dev_acc), "test set overlaps development data"


if __name__ == "__main__":
    main()
