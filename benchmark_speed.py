#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Wall-clock speed and scaling: pLIN vs MOB-suite (mob_typer) vs pling.

Each tool types random subsets (seed 7) of the 2,747 comparator-benchmark
plasmids, from FASTA files to labels, on the same machine.

  pLIN v4    plin_v4.py query (pyrodigal + exact/MMseqs2 family assignment +
             founder tree), --threads THREADS; timed with and without exact
             protein matching (--v4).
  pLIN       (v3) query mode, as in the app: load the released founder tree
             (data/plin_founder_tree_reference.npz, timed separately), then
             read each FASTA, compute its 4-mer vector and assign it on a
             session copy of the tree. Single-threaded (pLIN does not parallelise).
  MOB-suite  mob_typer per plasmid, THREADS plasmids in parallel.
  pling      pling cluster align --sourmash --cores THREADS, no visualisation.

Larger sizes are extrapolated linearly for pLIN and MOB-suite (constant
per-plasmid cost) and flagged as such. pling is not extrapolated: its cost
depends on how many plasmid pairs are related, not only on their number.

Timings are only meaningful on an otherwise idle machine, so the script stops
if another pling or mob_typer process is running (override with --force).
Run directories are never reused: an existing one aborts the script.

Usage:
  python benchmark_speed.py [--force] [--sizes 100 500 1000] [--pling-sizes 100 250 500]
                            [--inputs LIST --outdir DIR --v4]
Output: output/comparator_benchmark/speed/speed_results.tsv
"""

import argparse
import os
import random
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from Bio import SeqIO

from plin_founder import FounderTree, kmer_vector, session_copy

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CB = os.path.join(BASE_DIR, "output", "comparator_benchmark")
SP = os.path.join(CB, "speed")
TREE = os.path.join(BASE_DIR, "data", "plin_founder_tree_reference.npz")
MOB = os.path.expanduser("~/miniconda3/envs/mob_suite_env/bin/mob_typer")
PLING_BIN = os.path.expanduser("~/miniconda3/envs/pling_env/bin")
THREADS = 8
EXTRAPOLATE_TO = [10_000, 100_000]


def subset(paths, n):
    return sorted(random.Random(7).sample(paths, n))


def fresh_dir(path):
    if os.path.exists(path):
        sys.exit(f"{path} already exists; move it away before re-timing")
    os.makedirs(path)
    return path


def time_plin(paths, sizes):
    t = time.perf_counter()
    tree = FounderTree.from_npz(TREE)
    load_s = time.perf_counter() - t
    rows = [{"tool": "pLIN", "step": "load released tree (once per session)", "n_plasmids": 0,
             "threads": 1, "wall_s": round(load_s, 2), "source": "measured"}]
    for n in sizes:
        work = session_copy(tree)
        t = time.perf_counter()
        for p in subset(paths, n):
            work.assign(kmer_vector(str(next(SeqIO.parse(p, "fasta")).seq)), key=os.path.basename(p))
        rows.append({"tool": "pLIN", "step": "type plasmids", "n_plasmids": n, "threads": 1,
                     "wall_s": round(time.perf_counter() - t, 2), "source": "measured"})
        print(rows[-1], flush=True)
    return rows


def time_v4(paths, sizes, search_only):
    """pLIN v4 query mode as a separate process (includes loading the release, ~4 s).

    Evaluation plasmids are in the release, so whole-plasmid matching is
    always switched off; search_only also switches off exact protein matching
    (worst case for a plasmid unlike anything in the database)."""
    rows = []
    label = "pLIN v4 (all proteins searched)" if search_only else "pLIN v4 (exact protein matches used)"
    for n in sizes:
        d = fresh_dir(os.path.join(SP, f"v4_{'search' if search_only else 'typical'}_N{n}"))
        lst = os.path.join(d, "inputs.txt")
        open(lst, "w").write("\n".join(subset(paths, n)) + "\n")
        cmd = [sys.executable, os.path.join(BASE_DIR, "plin_v4.py"), "query", lst, "--out", os.path.join(d, "codes.tsv"),
               "--workdir", os.path.join(d, "work"), "--threads", str(THREADS), "--no-exact"]
        t = time.perf_counter()
        subprocess.run(cmd + (["--search-only"] if search_only else []), check=True, capture_output=True)
        rows.append({"tool": label, "step": "type plasmids", "n_plasmids": n, "threads": THREADS,
                     "wall_s": round(time.perf_counter() - t, 2), "source": "measured"})
        print(rows[-1], flush=True)
    return rows


def time_mob(paths, sizes):
    rows = []
    for n in sizes:
        out = fresh_dir(os.path.join(SP, f"mob_N{n}"))

        def run(p):
            o = os.path.join(out, os.path.basename(p)[:-len(".fasta")] + ".txt")
            subprocess.run([MOB, "-i", p, "-o", o, "-n", "1"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return os.path.exists(o) and os.path.getsize(o) > 0

        t = time.perf_counter()
        with ThreadPoolExecutor(THREADS) as ex:
            ok = sum(ex.map(run, subset(paths, n)))
        rows.append({"tool": "MOB-suite", "step": "type plasmids", "n_plasmids": n, "threads": THREADS,
                     "wall_s": round(time.perf_counter() - t, 2), "source": "measured", "failed": n - ok})
        print(rows[-1], flush=True)
    return rows


def time_pling(paths, sizes):
    rows = []
    env = dict(os.environ, PATH=PLING_BIN + os.pathsep + os.environ["PATH"])
    for n in sizes:
        d = fresh_dir(os.path.join(SP, f"pling_N{n}"))
        lst = os.path.join(d, "inputs.txt")
        open(lst, "w").write("\n".join(subset(paths, n)) + "\n")
        t = time.perf_counter()
        r = subprocess.run(["pling", "cluster", "align", "inputs.txt", "out", "--sourmash", "--cores", str(THREADS),
                            "--visualisation", "none"], cwd=d, env=env,
                           stdout=open(os.path.join(d, "pling.log"), "w"), stderr=subprocess.STDOUT)
        wall = time.perf_counter() - t
        if r.returncode != 0:
            sys.exit(f"pling failed for N={n}; see {d}/pling.log")
        rows.append({"tool": "pling", "step": "type plasmids", "n_plasmids": n, "threads": THREADS,
                     "wall_s": round(wall, 2), "source": "measured"})
        print(rows[-1], flush=True)
    return rows


def extrapolate(rows):
    out = []
    df = pd.DataFrame([r for r in rows if r["step"] == "type plasmids"])
    for tool, g in df.groupby("tool"):
        x, y = g.n_plasmids.values.astype(float), g.wall_s.values.astype(float)
        if tool == "pling":
            # Not extrapolated: pling's cost grows with the number of related pairs, so small random
            # subsets are dominated by fixed start-up time and a fitted power law is meaningless.
            continue
        else:                                      # per-plasmid cost is constant (pLIN v3/v4, MOB-suite)
            rate = y.sum() / x.sum()
            f, how = (lambda n: rate * n), "extrapolated (linear)"
        for n in EXTRAPOLATE_TO:
            out.append({"tool": tool, "step": "type plasmids", "n_plasmids": n, "threads": int(g.threads.iloc[0]),
                        "wall_s": round(float(f(n)), 1), "source": how})
    return out


def main():
    global SP
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--sizes", type=int, nargs="+", default=[100, 500, 1000])
    ap.add_argument("--pling-sizes", type=int, nargs="+", default=[100, 250, 500])
    ap.add_argument("--inputs", default=os.path.join(CB, "pling_inputs.txt"), help="FASTA path list to sample from")
    ap.add_argument("--outdir", default=SP)
    ap.add_argument("--v4", action="store_true", help="also time pLIN v4 query mode")
    args = ap.parse_args()
    SP = args.outdir
    busy = subprocess.run(["pgrep", "-f", "pling (cluster|add)|mob_typer"], capture_output=True).stdout
    if busy and not args.force:
        sys.exit("pling or mob_typer is running; timings would be distorted (use --force to run anyway)")

    os.makedirs(SP, exist_ok=True)
    paths = [l.strip() for l in open(args.inputs) if l.strip()]
    rows = time_plin(paths, args.sizes + [len(paths)])
    if args.v4:
        rows += time_v4(paths, args.sizes, search_only=False)
        rows += time_v4(paths, args.sizes, search_only=True)
    rows += time_mob(paths, args.sizes)
    rows += time_pling(paths, args.pling_sizes)
    rows += extrapolate(rows)
    res = pd.DataFrame(rows)
    res["s_per_plasmid"] = (res.wall_s / res.n_plasmids.replace(0, np.nan)).round(4)
    res.to_csv(os.path.join(SP, "speed_results.tsv"), sep="\t", index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
