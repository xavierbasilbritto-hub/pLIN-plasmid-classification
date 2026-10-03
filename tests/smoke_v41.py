#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Smoke test of pLIN v4.1 typing on any platform (run by the desktop build workflow on Windows, macOS
and Linux, and locally).

  1. MMseqs2 (bundled or installed) runs a protein search
  2. the v4.1 database release downloads and every file matches its published SHA-256 checksum
  3. two database plasmids typed from their sequence get their published codes

Usage:
  python tests/smoke_v41.py --release-dir DIR [--skip-download]
Exit code 0 when every check passes.
"""

import argparse
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from plin_v41_typer import PlinV41Release, download_release, find_mmseqs  # noqa: E402

EXPECTED = {"AB000510.1": "1.1.1.1.1.1", "AB011548.2": "2.2.2.2.2.2"}
PROT = ("MKKLLIAGLLAGSLLTGCATQPKHEEISAEAQRQIEQLKAELDALKAQNQQLRQELEDLKKALES"
        "AGFDVKLDETGRVLLTLPEDLLFASGSAELNSEGQAQLDALAAQLKQVDHGAHVVVQGHTDS")


def check_mmseqs(work):
    mm = find_mmseqs()
    assert mm, "MMseqs2 not found"
    print("MMseqs2:", mm)
    q, d = os.path.join(work, "q.faa"), os.path.join(work, "d.faa")
    open(q, "w").write(f">q1\n{PROT}\n")
    open(d, "w").write(f">t1\n{PROT}\n>t2\n{PROT[::-1]}\n")
    out = os.path.join(work, "hits.m8")
    r = subprocess.run([mm, "easy-search", q, d, out, os.path.join(work, "tmp"), "--format-output", "query,target,pident"],
                       capture_output=True, text=True)
    assert r.returncode == 0, f"mmseqs easy-search failed: {r.stderr[-500:]}"
    hits = open(out).read().split()
    assert hits[:2] == ["q1", "t1"], f"unexpected mmseqs output: {hits}"
    print("  ok  MMseqs2 search")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release-dir", required=True)
    ap.add_argument("--skip-download", action="store_true")
    args = ap.parse_args()
    work = tempfile.mkdtemp(prefix="plin_smoke_")
    check_mmseqs(work)
    if not args.skip_download:
        download_release(args.release_dir, progress=None)
        print("  ok  release downloaded and checksums verified")
    rel = PlinV41Release(args.release_dir, index_dir=os.path.join(work, "index"))
    from Bio import SeqIO
    recs = [(r.id, str(r.seq)) for r in SeqIO.parse(os.path.join(HERE, "data", "known_plasmids.fasta"), "fasta")]
    typed = rel.type_sequences(recs, threads=2).set_index("plasmid_id").pLIN_v41
    for acc, code in EXPECTED.items():
        assert typed[acc] == code, f"{acc}: got {typed[acc]}, expected {code}"
        print(f"  ok  {acc} -> {code}")
    print("ALL SMOKE TESTS PASSED")


if __name__ == "__main__":
    main()
