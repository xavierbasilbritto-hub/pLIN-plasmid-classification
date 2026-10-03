#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Consistency check for every document published with pLIN (run before each release).

Fails (exit code 1) when a published text file
  - quotes a number or claim from an earlier version (STALE below),
  - contains an em dash or an AI-tool marker,
  - or when a headline number in README.md differs from docs/PLIN_FACTS.json
    (generated from the result files by build_facts.py).
Pre-registration documents are frozen and skipped.

Usage:
  python check_docs.py
"""

import json
import os
import re
import subprocess
import sys

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
EXT = (".md", ".py", ".txt", ".cff", ".yml", ".yaml", ".sh", ".bat", ".command", ".spec", ".tsv")
SKIP = re.compile(r"(^|/)PREREGISTRATION|^output/|^reference/|^plasmid_sequences_for_training/|^for_review|^db_expansion")
STALE = {
    "133,326": "old database size", "79,305": "old database size", "72,959": "old database size",
    "72,556": "old database size", "82,477": "old code count", "57,886": "old code count",
    "1.1.2.4.7.13": "old single-linkage Swiss VIM-1 code", "ANI Equivalent": "levels are not ANI bands",
    "% ANI)": "levels are not ANI bands", "based on tetranucleotide (4-mer) composition distances": "v3 method",
    "Inc Groups Supported (20)": "28 groups",
}
ALLOW = {"plin_app.py": {"133,305"}, "check_docs.py": set(STALE) | {"—"}}
AI = re.compile(r"Co-Authored-By: Claude|Generated with \[Claude|noreply@anthropic|\U0001F916|\bChatGPT\b")
README_FACTS = ["db_plasmids", "clusters_L6", "lineage_F1", "lineage_F1_MOB", "lineage_F1_CI", "backbone_F1",
                "backbone_F1_MOB", "stability_runs", "requery_n", "speed_related_s", "speed_divergent_s",
                "speed_MOB_s", "replicon_cv_accuracy", "replicon_macro_F1", "test_pairs"]


def published_files():
    tracked = subprocess.run(["git", "ls-files"], cwd=BASE_DIR, capture_output=True, text=True).stdout.split()
    new = subprocess.run(["git", "ls-files", "--others", "--exclude-standard"], cwd=BASE_DIR,
                         capture_output=True, text=True).stdout.split()
    ignored = set(subprocess.run(["git", "ls-files", "--cached", "--ignored", "--exclude-standard"], cwd=BASE_DIR,
                                 capture_output=True, text=True).stdout.split())   # tracked but now excluded
    return sorted(f for f in (set(tracked) - ignored) | set(new) if f.endswith(EXT) and not SKIP.search(f)
                  and os.path.isfile(os.path.join(BASE_DIR, f)))


def main():
    problems = []
    for f in published_files():
        try:
            text = open(os.path.join(BASE_DIR, f), encoding="utf-8").read()
        except UnicodeDecodeError:
            continue
        allow = ALLOW.get(f, set())
        for pat, why in STALE.items():
            if pat in text and pat not in allow:
                problems.append(f"{f}: '{pat}' ({why})")
        if "—" in text and "—" not in allow:
            problems.append(f"{f}: {text.count(chr(0x2014))} em dash(es)")
        if f != "check_docs.py" and any(e in text for e in ("&mdash;", "&#8212;", "\\u2014")):
            problems.append(f"{f}: escaped em dash")
        if AI.search(text) and f != "check_docs.py":
            problems.append(f"{f}: AI-tool marker")
    facts = json.load(open(os.path.join(BASE_DIR, "docs", "PLIN_FACTS.json")))
    readme = open(os.path.join(BASE_DIR, "README.md"), encoding="utf-8").read()
    for k in README_FACTS:
        if facts[k]["text"] not in readme:
            problems.append(f"README.md: fact {k} = '{facts[k]['text']}' not found")
    if problems:
        print("\n".join(problems))
        print(f"\nFAILED: {len(problems)} problem(s)")
        sys.exit(1)
    print(f"OK: {len(published_files())} published files checked, README matches {len(README_FACTS)} facts")


if __name__ == "__main__":
    main()
