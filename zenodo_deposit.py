#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
Deposit a pLIN database release on Zenodo and reserve its DOI.

Steps:
  1. check every release file against the SHA-256 checksums in DATABASE_VERSION.json
  2. create a Zenodo deposition with the metadata below and a reserved DOI
  3. upload the files and check that Zenodo's MD5 for each file matches the local one
  4. stop as a draft (nothing is public); with --publish, publish an existing draft

The Zenodo personal access token (scopes deposit:write and deposit:actions) is read from
the ZENODO_TOKEN environment variable or the file ~/.zenodo_token; it is never printed.

Usage:
  python zenodo_deposit.py [--release-dir DIR] [--sandbox]      create the draft, reserve the DOI
  python zenodo_deposit.py --publish [--sandbox]                 publish the draft (permanent)
  python zenodo_deposit.py --code v4.1.2 [--publish-code]        archive the source code of a git tag
Output: zenodo_deposition.json (deposition id, reserved DOI, links); zenodo_code_<tag>.json for --code
"""

import argparse
import hashlib
import json
import os
import sys

import requests

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
STATE = os.path.join(BASE_DIR, "zenodo_deposition.json")
REPO = "https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification"

UMCG = ("Department of Medical Microbiology and Infection Prevention, University Medical Center Groningen, "
        "University of Groningen, Groningen, The Netherlands")
CREATORS = [
    {"name": "Xavier, Basil Britto", "affiliation": UMCG, "orcid": "0000-0002-5897-0240"},
    {"name": "Bari, Anurag Kumar", "affiliation": UMCG},
    {"name": "Sinha, Bhanu", "affiliation": UMCG},
    {"name": "Rossen, John W. A.", "affiliation": UMCG + "; Isala Hospital, Zwolle, The Netherlands; "
                                                       "University of Utah School of Medicine, Salt Lake City, UT, USA"},
]


def token():
    t = os.environ.get("ZENODO_TOKEN")
    if not t and os.path.exists(os.path.expanduser("~/.zenodo_token")):
        t = open(os.path.expanduser("~/.zenodo_token")).read().strip()
    if not t:
        sys.exit("No Zenodo token: set ZENODO_TOKEN or write it to ~/.zenodo_token (scopes deposit:write, deposit:actions)")
    return t


def digest(path, algo):
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def metadata(dv):
    lv = "".join(f"<li>{k}: {v}</li>" for k, v in dv["levels"].items())
    files = "".join(f"<li>{f} (SHA-256 {s})</li>" for f, s in dv["sha256"].items())
    desc = (
        f"<p>Database release {dv['database_version']} of pLIN (plasmid lineage identification number), a permanent "
        f"six-level nomenclature for bacterial plasmids. The release assigns pLIN {dv['scheme'].split()[-1]} codes to "
        f"{dv['unique_plasmids']:,} unique plasmid sequences from PLSDB and NCBI and contains everything needed to type new "
        f"plasmids with the pLIN software: the codes, the k-mer sketch index, the protein-family catalogue and the "
        f"sequence-hash table.</p><p>Levels:</p><ul>{lv}</ul>"
        f"<p>Database codes never change in later releases: new plasmids are added after the existing ones.</p>"
        f"<p>Files:</p><ul>{files}<li>DATABASE_VERSION.json (release record and checksums)</li></ul>"
        f"<p>Software, documentation and the pre-registered evaluation: <a href=\"{REPO}\">{REPO}</a>. "
        f"Place the files in data/plin_v41 (or ~/.plin/plin_v41) for the pLIN app and command-line typer.</p>")
    return {
        "title": f"pLIN {dv['scheme'].split()[-1]} plasmid nomenclature database, release {dv['database_version']}",
        "upload_type": "dataset",
        "description": desc,
        "creators": CREATORS,
        "access_right": "open",
        "license": "cc-by-4.0",
        "version": dv["database_version"],
        "keywords": ["plasmid", "plasmid typing", "nomenclature", "pLIN", "antimicrobial resistance",
                     "genomic epidemiology", "One Health"],
        "related_identifiers": [
            {"identifier": REPO, "relation": "isSupplementTo", "resource_type": "software"},
            {"identifier": "10.21203/rs.3.rs-10481391/v1", "relation": "isDescribedBy", "resource_type": "publication-preprint"},
        ],
        "prereserve_doi": True,
    }


def create(api, headers, release_dir):
    dv = json.load(open(os.path.join(release_dir, "DATABASE_VERSION.json")))
    files = list(dv["sha256"]) + ["DATABASE_VERSION.json"]
    print("Checking files against DATABASE_VERSION.json ...")
    md5 = {}
    for f in files:
        p = os.path.join(release_dir, f)
        if f in dv["sha256"] and digest(p, "sha256") != dv["sha256"][f]:
            sys.exit(f"SHA-256 mismatch for {f}; not uploading")
        md5[f] = digest(p, "md5")
        print(f"  ok  {f}")
    r = requests.post(f"{api}/deposit/depositions", headers=headers, json={"metadata": metadata(dv)}, timeout=60)
    r.raise_for_status()
    dep = r.json()
    doi = dep["metadata"]["prereserve_doi"]["doi"]
    json.dump({"id": dep["id"], "doi": doi, "links": dep["links"], "sandbox": "sandbox" in api, "published": False},
              open(STATE, "w"), indent=1)
    print(f"Draft deposition {dep['id']} created; reserved DOI {doi}")
    bucket = dep["links"]["bucket"]
    for f in files:
        size = os.path.getsize(os.path.join(release_dir, f))
        print(f"Uploading {f} ({size / 1e6:.0f} MB) ...")
        with open(os.path.join(release_dir, f), "rb") as fh:
            r = requests.put(f"{bucket}/{f}", data=fh, headers=headers, timeout=None)
        r.raise_for_status()
        got = r.json().get("checksum", "").replace("md5:", "")
        if got != md5[f]:
            sys.exit(f"Checksum mismatch after upload of {f} (Zenodo {got}, local {md5[f]}); the draft is not published")
        print(f"  ok  {f} (MD5 verified)")
    print(f"\nDraft complete and NOT public. Review it at {dep['links']['html']}")
    print(f"Reserved DOI: {doi}  (becomes active when the draft is published)")


def publish(api, headers):
    st = json.load(open(STATE))
    r = requests.post(f"{api}/deposit/depositions/{st['id']}/actions/publish", headers=headers, timeout=120)
    r.raise_for_status()
    st["published"] = True
    st["record"] = r.json()["links"].get("record_html") or r.json()["links"].get("html")
    json.dump(st, open(STATE, "w"), indent=1)
    print(f"Published. DOI {st['doi']}: https://doi.org/{st['doi']}")


def code_deposit(api, headers, tag, publish_now):
    """Archive the source code of a git tag (git archive) as a Zenodo software record."""
    import subprocess
    import tempfile
    version = tag.lstrip("v")
    name = f"pLIN-plasmid-classification-{version}.zip"
    path = os.path.join(tempfile.mkdtemp(), name)
    subprocess.run(["git", "archive", "--format=zip", f"--prefix=pLIN-plasmid-classification-{version}/", "-o", path, tag],
                   cwd=BASE_DIR, check=True)
    commit = subprocess.run(["git", "rev-list", "-n", "1", tag], cwd=BASE_DIR, capture_output=True, text=True).stdout.strip()
    meta = {
        "title": f"pLIN: permanent plasmid nomenclature software, version {version}",
        "upload_type": "software",
        "description": (f"<p>Source code of pLIN {version} (plasmid lineage identification number), git tag {tag} "
                        f"(commit {commit}) of <a href=\"{REPO}\">{REPO}</a>: the typing software, desktop application, "
                        f"database build, and the scripts of the pre-registered confirmatory evaluation, the registered "
                        f"external validation and the sensitivity analyses.</p><p>The matching database release "
                        f"db-2026.10.03 is archived at https://doi.org/10.5281/zenodo.23126057.</p>"),
        "creators": CREATORS, "access_right": "open", "license": "gpl-3.0-or-later", "version": version,
        "keywords": ["plasmid", "plasmid typing", "nomenclature", "pLIN", "antimicrobial resistance", "software"],
        "related_identifiers": [
            {"identifier": f"{REPO}/tree/{tag}", "relation": "isIdenticalTo", "resource_type": "software"},
            {"identifier": "10.5281/zenodo.23126057", "relation": "isSupplementedBy", "resource_type": "dataset"}],
        "prereserve_doi": True,
    }
    r = requests.post(f"{api}/deposit/depositions", headers=headers, json={"metadata": meta}, timeout=60)
    r.raise_for_status()
    dep = r.json()
    doi = dep["metadata"]["prereserve_doi"]["doi"]
    with open(path, "rb") as fh:
        u = requests.put(f"{dep['links']['bucket']}/{name}", data=fh, headers=headers, timeout=None)
    u.raise_for_status()
    if u.json().get("checksum", "").replace("md5:", "") != digest(path, "md5"):
        sys.exit("Checksum mismatch after upload; the draft is not published")
    state = {"id": dep["id"], "doi": doi, "tag": tag, "commit": commit, "published": False}
    if publish_now:
        p = requests.post(f"{api}/deposit/depositions/{dep['id']}/actions/publish", headers=headers, timeout=120)
        p.raise_for_status()
        state["published"] = True
    json.dump(state, open(os.path.join(BASE_DIR, f"zenodo_code_{tag}.json"), "w"), indent=1)
    print(f"Code {tag} ({commit[:7]}) {'published' if publish_now else 'drafted'}: DOI {doi}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release-dir", default=os.path.join(BASE_DIR, "output", "backbone_v41", "release"))
    ap.add_argument("--sandbox", action="store_true", help="use sandbox.zenodo.org (test DOIs)")
    ap.add_argument("--publish", action="store_true", help="publish the draft recorded in zenodo_deposition.json")
    ap.add_argument("--code", metavar="TAG", help="archive the source code of this git tag instead of the database")
    ap.add_argument("--publish-code", action="store_true", help="with --code: publish immediately (permanent)")
    args = ap.parse_args()
    api = "https://sandbox.zenodo.org/api" if args.sandbox else "https://zenodo.org/api"
    headers = {"Authorization": f"Bearer {token()}"}
    if args.code:
        code_deposit(api, headers, args.code, args.publish_code)
        return
    if args.publish:
        publish(api, headers)
    else:
        if os.path.exists(STATE) and not json.load(open(STATE)).get("published"):
            sys.exit(f"A draft already exists ({STATE}); publish it with --publish or delete it on Zenodo first")
        create(api, headers, args.release_dir)


if __name__ == "__main__":
    main()
