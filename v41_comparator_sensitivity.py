#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Comparator parameter sensitivity on the confirmatory test plasmids (secondary analysis, not pre-registered).

The registered comparison ran pling and mge-cluster with default settings. Here each is rerun on the same
2,000 test plasmids with other documented settings, and scored on the same 3,337 weighted pairs, so that
the comparison with pLIN uses each tool's best setting:
  pling        reuse of the registered run's distances (pling cluster skip --reuse_previous) with
               DCJ-Indel thresholds 1, 2, 3, 4 (default), 6, 8, 10 at containment 0.5 (default), and
               containment distance 0.3 and 0.4 at DCJ 4 (pling requires the reused run's containment
               threshold to be at least the new one)
  mge-cluster  new models on the 2,000 plasmids with perplexity 10, 30 (default), 50 and minimum cluster
               size 5, 10, 30 (default)

Usage:
  python v41_comparator_sensitivity.py
Output: output/backbone_v41/confirm/comparator_sensitivity.tsv
"""

import glob
import os
import shutil
import subprocess

import numpy as np
import pandas as pd

from v41_confirm_eval import C
from v4_evaluate import TRUTHS, metrics, same_group

PLING = os.path.expanduser("~/miniconda3/envs/pling_env/bin/pling")
MGE = os.path.expanduser("~/miniconda3/envs/mge_cluster_env/bin/mge_cluster")
ENV = dict(os.environ, PATH=os.path.dirname(PLING) + os.pathsep + os.environ["PATH"])


def pling_labels(out):
    t = glob.glob(os.path.join(out, "dcj_thresh_*_graph", "objects", "typing.tsv"))
    c = os.path.join(out, "containment", "containment_communities", "objects", "communities.tsv")
    if not t or not os.path.exists(c):
        return None
    sub = dict(pd.read_csv(t[0], sep="\t", dtype=str).iloc[:, :2].values)
    com = dict(pd.read_csv(c, sep="\t", dtype=str).iloc[:, :2].values)
    return sub, com


def run_pling(dcj, cont):
    out = os.path.join(C, "pling_sens", f"dcj{dcj}_cont{cont}")
    if pling_labels(out) is None:
        prev = os.path.join(C, "pling", "out")
        # pling 3.0.2 --reuse_previous copies the previous distances into <out>/containment and reads the
        # alignment batches from <out>/batches, but creates neither: prepare both
        shutil.rmtree(out, ignore_errors=True)
        os.makedirs(os.path.join(out, "containment"))
        shutil.copytree(os.path.join(prev, "batches"), os.path.join(out, "batches"))
        subprocess.run([PLING, "cluster", "skip", os.path.join(C, "pling", "inputs.txt"), out, prev,
                        "--reuse_previous", prev, "--dcj", str(dcj), "--containment_distance", str(cont),
                        "--visualisation", "none", "--cores", "8"], env=ENV, capture_output=True, cwd=os.path.join(C, "pling"))
    return pling_labels(out)


def run_mge(perp, minc):
    out = os.path.join(C, "mge_sens", f"perp{perp}_min{minc}")
    res = os.path.join(out, "mge-cluster_results.csv")
    if not os.path.exists(res):
        # the unitig table depends only on the input plasmids: reuse it (mge-cluster skips unitig calling when
        # <outdir>/mge-cluster.rtab exists), so that only the t-SNE embedding and HDBSCAN are rerun
        shared = glob.glob(os.path.join(C, "mge_sens", "*", "sorted_mge-cluster.rtab"))
        os.makedirs(out, exist_ok=True)
        if shared and not os.path.exists(os.path.join(out, "sorted_mge-cluster.rtab")):
            src = os.path.dirname(shared[0])
            for f in ("mge-cluster.rtab", "sorted_mge-cluster.rtab"):
                os.symlink(os.path.join(src, f), os.path.join(out, f))
        subprocess.run([MGE, "--create", "--input", os.path.join(C, "mge", "inputs.txt"), "--outdir", out,
                        "--perplexity", str(perp), "--min_cluster", str(minc), "--threads", "8"], capture_output=True)
    if not os.path.exists(res):
        return None
    m = pd.read_csv(res, dtype=str)
    return {os.path.basename(s).rsplit(".", 1)[0] if s.endswith((".fasta", ".fa", ".fna")) else s:
            (None if c.strip() == "-1" else c.strip()) for s, c in zip(m.Sample_Name, m.Standard_Cluster)}


def main():
    pairs = pd.read_csv(os.path.join(C, "truth_pairs.tsv"), sep="\t").reset_index(drop=True)
    a, b, w = pairs.plasmid_shorter, pairs.plasmid_longer, pairs.weight.values
    rows = []

    ids = set(a) | set(b)

    def score(tool, setting, lab):
        lab = {i: lab.get(i) for i in ids}                # plasmids a tool leaves out are unassigned
        pred = same_group(lab, a, b)
        for t in TRUTHS:
            p, r, f = metrics(pred, pairs[t].values, w)
            rows.append({"tool": tool, "setting": setting, "truth": t, "precision_w": p, "recall_w": r, "F1_w": f,
                         "assigned": float(np.mean([lab.get(i) is not None for i in ids]))})

    for dcj, cont in [(1, 0.5), (2, 0.5), (3, 0.5), (4, 0.5), (6, 0.5), (8, 0.5), (10, 0.5), (4, 0.4), (4, 0.3)]:
        lab = run_pling(dcj, cont)
        if lab is None:
            print("pling failed", dcj, cont)
            continue
        default = " (default)" if (dcj, cont) == (4, 0.5) else ""
        score("pling subcommunity", f"dcj {dcj}, containment {cont}{default}", lab[0])
        score("pling community", f"dcj {dcj}, containment {cont}{default}", lab[1])
        print("pling done", dcj, cont, flush=True)
    for perp in (10, 30, 50):
        for minc in (5, 10, 30):
            lab = run_mge(perp, minc)
            if lab is None:
                print("mge-cluster failed", perp, minc)
                continue
            default = " (default)" if (perp, minc) == (30, 30) else ""
            score("mge-cluster", f"perplexity {perp}, min cluster {minc}{default}", lab)
            print("mge-cluster done", perp, minc, flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(C, "comparator_sensitivity.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(res.pivot_table(index=["tool", "setting"], columns="truth", values="F1_w").round(3).to_string())


if __name__ == "__main__":
    main()
