#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
NCBI record-creation date for every release plasmid and every outbreak plasmid.

PLSDB 2024_05_31_v2 metadata (NUCCORE_CreateDate) where available; otherwise the
NCBI nuccore "CreateDate" from E-utilities esummary (batches of 200, cached).
These are deposition dates, not sample-collection dates.

Usage:
  python v41_dates.py
Output: output/backbone_v41/accession_dates.tsv (accession, create_date, source)
"""

import json
import os
import time
import urllib.parse
import urllib.request

import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "output", "backbone_v41", "accession_dates.tsv")
CACHE = os.path.join(BASE_DIR, "output", "backbone_v41", "ncbi_esummary_cache.json")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def esummary(accs):
    data = urllib.parse.urlencode({"db": "nuccore", "id": ",".join(accs), "retmode": "json"}).encode()
    for attempt in range(5):
        try:
            r = json.loads(urllib.request.urlopen(EUTILS, data=data, timeout=120).read())
            res = r.get("result", {})
            out = {}
            for uid in res.get("uids", []):
                d = res[uid]
                acc = d.get("accessionversion", "")
                out[acc] = out[acc.split(".")[0]] = d.get("createdate", "")   # versioned and bare accession
            return out
        except Exception:
            time.sleep(2 * (attempt + 1))
    return {}


def main():
    rel = pd.read_csv(os.path.join(BASE_DIR, "output", "backbone_v41", "release", "plin_v41_codes.tsv.gz"),
                      sep="\t", usecols=["accession"])
    ob = pd.read_csv(os.path.join(BASE_DIR, "output", "outbreak_validation_founder_results.tsv"), sep="\t")
    accs = sorted(set(rel.accession) | set(ob.accession))
    plsdb = pd.read_csv(os.path.join(BASE_DIR, "output", "backbone_v4", "onehealth", "plsdb_meta", "nuccore.csv"),
                        usecols=["NUCCORE_ACC", "NUCCORE_CreateDate"])
    dates = {a: (d, "PLSDB") for a, d in zip(plsdb.NUCCORE_ACC, plsdb.NUCCORE_CreateDate)}
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    todo = [a for a in accs if a not in dates and a not in cache]
    for k in range(0, len(todo), 200):
        got = esummary(todo[k:k + 200])
        for a in todo[k:k + 200]:
            cache[a] = got.get(a, "")
        if (k // 200) % 25 == 0:
            json.dump(cache, open(CACHE, "w"))
            print(f"  {k + 200:,}/{len(todo):,}", flush=True)
        time.sleep(0.35)
    json.dump(cache, open(CACHE, "w"))
    rows = []
    for a in accs:
        if a in dates:
            d, src = dates[a]
        else:
            d, src = cache.get(a, ""), "NCBI esummary"
        rows.append({"accession": a, "create_date": pd.to_datetime(d, errors="coerce"), "source": src})
    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep="\t", index=False)
    print(f"{len(df):,} accessions; dated {df.create_date.notna().sum():,}; missing {df.create_date.isna().sum():,}")


if __name__ == "__main__":
    main()
