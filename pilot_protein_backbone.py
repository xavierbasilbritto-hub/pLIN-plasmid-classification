#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Pilot: does backbone protein-family sharing add information beyond 4-mer
composition for judging plasmid relatedness?

Uses the plasmid pairs already aligned by validate_alignment_backbone.py
(single-linkage and founder runs pooled) as the truth set.

  1. Proteins predicted with Prodigal (-p meta) for every plasmid in the pairs.
  2. Clustered with MMseqs2 easy-cluster (>= 50% identity, >= 80% coverage)
     into protein families.
  3. Proteins overlapping annotated cargo (AMRFinderPlus AMR elements, curated
     IS hits, composite-transposon spans) are cargo; the rest are backbone.
  4. Per pair: backbone family containment (share of the shorter plasmid's
     backbone families present in the longer) and Jaccard.
  5. Evaluation against alignment truth (same lineage: bidirectional AF >= 0.8
     and identity >= 99%; related backbone: shorter-plasmid backbone AF >= 0.5):
     ROC AUC of cosine distance vs protein containment, and whether protein
     sharing recovers related pairs that the founder L6 code misses.

Usage:
  python pilot_protein_backbone.py --threads 12
Output: output/protein_pilot/{pair_protein_metrics.tsv, summary.json}
"""

import os
import json
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from validate_alignment_backbone import fasta_paths, load_accessory_intervals

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "output", "protein_pilot")
MMSEQS = os.path.expanduser("~/miniconda3/envs/pLIN_tools/bin/mmseqs")
PAIRS = [os.path.join(BASE_DIR, "output", "alignment_validation_single_linkage_v3.2.3", "pair_metrics.tsv"),
         os.path.join(BASE_DIR, "output", "alignment_validation", "pair_metrics.tsv")]
FOUNDER = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")


def prodigal(args):
    pid, path = args
    faa = os.path.join(OUT, "faa", f"{pid}.faa")
    if not os.path.exists(faa):
        subprocess.run(["prodigal", "-i", path, "-a", faa, "-p", "meta", "-q", "-o", os.devnull], check=True)
    return pid, faa


def read_faa(pid, faa):
    """Yield (protein_id, start, end) from Prodigal headers: >id # start # end # strand # ..."""
    for line in open(faa):
        if line.startswith(">"):
            parts = line[1:].split(" # ")
            yield parts[0].split()[0], int(parts[1]), int(parts[2])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=12)
    args = ap.parse_args()
    os.makedirs(os.path.join(OUT, "faa"), exist_ok=True)

    pairs = pd.concat([pd.read_csv(p, sep="\t") for p in PAIRS]).drop_duplicates(["plasmid_shorter", "plasmid_longer"])
    ids = sorted(set(pairs.plasmid_shorter) | set(pairs.plasmid_longer))
    paths = fasta_paths()
    print(f"[1/4] Prodigal on {len(ids)} plasmids ...", flush=True)
    with ThreadPoolExecutor(args.threads) as ex:
        faas = dict(ex.map(prodigal, [(i, paths[i]) for i in ids]))

    # combined fasta with plasmid-tagged ids; record cargo/backbone status
    intervals, _ = load_accessory_intervals()
    combined = os.path.join(OUT, "all_proteins.faa")
    prot_owner, prot_cargo = {}, {}
    with open(combined, "w") as out:
        for pid in ids:
            cargo_iv = intervals.get(pid, [])
            coords = {p: (s, e) for p, s, e in read_faa(pid, faas[pid])}
            for line in open(faas[pid]):
                if line.startswith(">"):
                    raw = line[1:].split()[0]
                    tag = f"{pid}|{raw}"
                    s, e = coords[raw]
                    overlap = any(min(e, b) - max(s, a) + 1 >= 0.5 * (e - s + 1) for a, b in cargo_iv)
                    prot_owner[tag], prot_cargo[tag] = pid, overlap
                    out.write(f">{tag}\n")
                else:
                    out.write(line)
    print(f"  {len(prot_owner):,} proteins, {sum(prot_cargo.values()):,} cargo", flush=True)

    print("[2/4] MMseqs2 clustering (50% id, 80% cov) ...", flush=True)
    prefix = os.path.join(OUT, "clu")
    if not os.path.exists(prefix + "_cluster.tsv"):
        subprocess.run([MMSEQS, "easy-cluster", combined, prefix, os.path.join(OUT, "tmp"),
                        "--min-seq-id", "0.5", "-c", "0.8", "--cov-mode", "0", "--threads", str(args.threads)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    fam = dict(pd.read_csv(prefix + "_cluster.tsv", sep="\t", header=None, names=["rep", "member"])
               .set_index("member")["rep"])
    backbone_fams, all_fams = {}, {}
    for tag, pid in prot_owner.items():
        f = fam[tag]
        all_fams.setdefault(pid, set()).add(f)
        if not prot_cargo[tag]:
            backbone_fams.setdefault(pid, set()).add(f)
    print(f"  {len(set(fam.values())):,} protein families", flush=True)

    print("[3/4] Pair metrics ...", flush=True)
    founder = pd.read_csv(FOUNDER, sep="\t").drop_duplicates("plasmid_id").set_index("plasmid_id").pLIN

    def shared(a, b):
        n = 0
        for x, y in zip(a.split("."), b.split(".")):
            if x != y:
                break
            n += 1
        return n

    rows = []
    for _, r in pairs.iterrows():
        a, b = r.plasmid_shorter, r.plasmid_longer
        ba, bb = backbone_fams.get(a, set()), backbone_fams.get(b, set())
        fa, fb = all_fams.get(a, set()), all_fams.get(b, set())
        rows.append({
            "plasmid_shorter": a, "plasmid_longer": b, "cosine_distance": r.cosine_distance,
            "AF_min": r.AF_min, "AF_shorter": r.AF_shorter, "ANI_aln": r.ANI_aln, "bb_AF_shorter": r.bb_AF_shorter,
            "n_bb_fam_shorter": len(ba),
            "bb_containment": len(ba & bb) / len(ba) if ba else np.nan,
            "bb_jaccard": len(ba & bb) / len(ba | bb) if (ba | bb) else np.nan,
            "all_containment": len(fa & fb) / len(fa) if fa else np.nan,
            "founder_shared_levels": shared(founder[a], founder[b]),
        })
    m = pd.DataFrame(rows)
    m["same_lineage"] = (m.AF_min >= 0.8) & (m.ANI_aln >= 99)
    m["related_backbone"] = m.bb_AF_shorter >= 0.5
    m.to_csv(os.path.join(OUT, "pair_protein_metrics.tsv"), sep="\t", index=False)

    print("[4/4] Evaluation ...", flush=True)
    s = {"n_pairs": int(len(m)), "n_plasmids": len(ids), "n_families": int(len(set(fam.values())))}
    ok = m.dropna(subset=["bb_containment", "ANI_aln"])
    for truth in ["same_lineage", "related_backbone"]:
        t = ok[truth]
        s[f"AUC_{truth}"] = {
            "cosine_distance (lower=related)": round(roc_auc_score(t, -ok.cosine_distance), 4),
            "backbone_protein_containment": round(roc_auc_score(t, ok.bb_containment), 4),
            "backbone_protein_jaccard": round(roc_auc_score(t, ok.bb_jaccard), 4),
            "all_protein_containment": round(roc_auc_score(t, ok.all_containment), 4),
        }
    # coarse levels: among pairs sharing only L1-L4 (founder), does protein sharing separate related backbones?
    coarse = ok[ok.founder_shared_levels.between(1, 4)]
    if coarse.related_backbone.nunique() == 2:
        s["coarse_L1_L4_pairs"] = {
            "n": int(len(coarse)), "n_related_backbone": int(coarse.related_backbone.sum()),
            "AUC_cosine": round(roc_auc_score(coarse.related_backbone, -coarse.cosine_distance), 4),
            "AUC_backbone_containment": round(roc_auc_score(coarse.related_backbone, coarse.bb_containment), 4)}
    # recovery of founder L6 misses
    l6 = ok.founder_shared_levels >= 6
    rel = ok.same_lineage
    for thr in (0.8, 0.9, 0.95):
        for cos_cap in (0.005, 0.01):
            extra = (~l6) & (ok.bb_containment >= thr) & (ok.cosine_distance <= cos_cap)
            pred = l6 | extra
            s[f"L6_or_protein>={thr}_cos<={cos_cap}"] = {
                "recall": round(float((pred & rel).sum() / rel.sum()), 4),
                "precision": round(float((pred & rel).sum() / pred.sum()), 4),
                "pairs_added": int(extra.sum()), "added_that_are_same_lineage": int((extra & rel).sum())}
    s["founder_L6_alone"] = {"recall": round(float((l6 & rel).sum() / rel.sum()), 4),
                             "precision": round(float((l6 & rel).sum() / l6.sum()), 4)}
    json.dump(s, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
