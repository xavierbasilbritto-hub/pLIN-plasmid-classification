#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
"""
Download Gram-positive plasmid sequences for pLIN classifier training.

Expands the Inc group classifier from 20 Gram-negative groups to include
S. aureus and E. faecium plasmid rep types, addressing limitations L4
(database bias) and L6 (limited to 20 Inc groups).

Rep type naming follows Lozano et al. 2012 (Staphylococcus) and
Jensen et al. 2010 / Clewell 2011 (Enterococcus).

Usage:
    python download_gram_positive_training.py
"""

import os
import sys
import time
import subprocess
import xml.etree.ElementTree as ET

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TRAINING_DIR = os.path.join(BASE_DIR, "plasmid_sequences_for_training")

# Entrez base URL
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

# ── Well-characterized reference plasmids per rep type ──────────────────────
# Each dict: rep_type -> {description, organism, search_query, known_accessions}
#
# known_accessions: manually curated from literature (gold standard)
# search_query: NCBI Nucleotide search to find additional examples
#
# S. aureus rep types (Lozano et al. 2012, Jensen et al. 2006)
# E. faecium rep types (Jensen et al. 2010, Palmer et al. 2010)

GRAM_POSITIVE_GROUPS = {
    # ── S. aureus plasmid families ──────────────────────────────────
    "rep5_pT181": {
        "description": "pT181-family small rolling-circle plasmids (tetracycline resistance)",
        "organism": "Staphylococcus aureus",
        "known_accessions": [
            "J01764.1",    # pT181 prototype
            "M19465.1",    # pC194
            "V01277.1",    # pUB110
            "M36227.1",    # pE194
            "U40259.1",    # pSK639
        ],
        "search_query": (
            '("Staphylococcus aureus"[Organism]) AND plasmid[Title] '
            'AND ("2000"[SLEN] : "6000"[SLEN]) AND refseq[Filter]'
        ),
        "max_download": 50,
    },
    "rep7_pSK1": {
        "description": "pSK1-family multi-resistance plasmids",
        "organism": "Staphylococcus aureus",
        "known_accessions": [
            "GQ900391.1",  # pSK1
            "AF203376.1",  # pSK41
            "GQ900381.1",  # pSK4
            "AF077865.1",  # pLW043
        ],
        "search_query": (
            '("Staphylococcus aureus"[Organism]) AND plasmid[Title] '
            'AND ("20000"[SLEN] : "50000"[SLEN]) AND refseq[Filter]'
        ),
        "max_download": 50,
    },
    "rep10_pI258": {
        "description": "pI258-family heavy metal resistance plasmids",
        "organism": "Staphylococcus aureus",
        "known_accessions": [
            "GQ900398.1",  # pI258
            "AP003139.1",  # pN315 (MRSA plasmid)
            "BA000018.3",  # pVRSA (vancomycin resistance)
        ],
        "search_query": (
            '("Staphylococcus aureus"[Organism]) AND plasmid[Title] '
            'AND ("25000"[SLEN] : "40000"[SLEN]) AND refseq[Filter]'
        ),
        "max_download": 50,
    },
    "rep16_SAP": {
        "description": "SAP-type (Staphylococcal Accessory Plasmid) small plasmids",
        "organism": "Staphylococcus aureus",
        "known_accessions": [
            "AY355285.1",  # SAP099B
            "GQ900404.1",  # pSAS
        ],
        "search_query": (
            '("Staphylococcus aureus"[Organism]) AND plasmid[Title] '
            'AND ("1000"[SLEN] : "5000"[SLEN]) AND refseq[Filter]'
        ),
        "max_download": 50,
    },
    "rep20_pSK41": {
        "description": "pSK41-family conjugative multiresistance plasmids",
        "organism": "Staphylococcus aureus",
        "known_accessions": [
            "AF203376.1",  # pSK41
            "AF411935.1",  # pGO1
            "GQ900406.1",  # pUSA03
        ],
        "search_query": (
            '("Staphylococcus aureus"[Organism]) AND plasmid[Title] '
            'AND ("40000"[SLEN] : "100000"[SLEN]) AND refseq[Filter]'
        ),
        "max_download": 50,
    },
    "rep19_pWBG": {
        "description": "pWBG749-family beta-lactamase plasmids",
        "organism": "Staphylococcus aureus",
        "known_accessions": [
            "GQ900405.1",  # pWBG749
            "GQ900393.1",  # pSK68
        ],
        "search_query": (
            '("Staphylococcus aureus"[Organism]) AND plasmid[Title] '
            'AND ("3000"[SLEN] : "15000"[SLEN]) AND refseq[Filter]'
        ),
        "max_download": 50,
    },

    # ── E. faecium plasmid families ────────────────────────────────
    "rep9_pAD1": {
        "description": "pAD1-family pheromone-responsive conjugative plasmids",
        "organism": "Enterococcus faecalis",
        "known_accessions": [
            "L01794.1",    # pAD1
            "AF394225.1",  # pCF10
            "X92945.2",    # pPD1
        ],
        "search_query": (
            '("Enterococcus"[Organism]) AND plasmid[Title] '
            'AND ("50000"[SLEN] : "80000"[SLEN]) AND refseq[Filter]'
        ),
        "max_download": 50,
    },
    "rep11_pCF10": {
        "description": "pCF10-family tetracycline conjugative plasmids",
        "organism": "Enterococcus faecalis",
        "known_accessions": [
            "AF394225.1",  # pCF10
            "AY855841.1",  # pTEF1
        ],
        "search_query": (
            '("Enterococcus faecalis"[Organism]) AND plasmid[Title] '
            'AND ("60000"[SLEN] : "80000"[SLEN]) AND refseq[Filter]'
        ),
        "max_download": 50,
    },
    "rep14_pRUM": {
        "description": "pRUM-family vancomycin resistance transfer plasmids",
        "organism": "Enterococcus faecium",
        "known_accessions": [
            "EF507804.1",  # pRUM
            "AF521699.1",  # pLG2
        ],
        "search_query": (
            '("Enterococcus faecium"[Organism]) AND plasmid[Title] '
            'AND refseq[Filter]'
        ),
        "max_download": 50,
    },
    "rep17_pRE25": {
        "description": "pRE25-family multi-resistance broad-host-range plasmids",
        "organism": "Enterococcus faecalis",
        "known_accessions": [
            "X92946.1",    # pRE25
            "AY234334.1",  # pHTbeta prototype
        ],
        "search_query": (
            '("Enterococcus"[Organism]) AND plasmid[Title] '
            'AND ("30000"[SLEN] : "60000"[SLEN]) AND refseq[Filter]'
        ),
        "max_download": 50,
    },
    "rep18_pHTbeta": {
        "description": "pHTbeta-family small mobilizable plasmids",
        "organism": "Enterococcus faecium",
        "known_accessions": [
            "AY234334.1",  # pHTbeta
        ],
        "search_query": (
            '("Enterococcus faecium"[Organism]) AND plasmid[Title] '
            'AND ("5000"[SLEN] : "30000"[SLEN]) AND refseq[Filter]'
        ),
        "max_download": 50,
    },
}


def _curl_fetch(url, timeout=30):
    """Fetch URL content using curl (more reliable than urllib in some envs)."""
    try:
        result = subprocess.run(
            ["curl", "-s", "-f", "--max-time", str(timeout), url],
            capture_output=True, text=True, timeout=timeout + 5,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout
    except Exception:
        pass
    return None


def fetch_fasta(accession, retries=3):
    """Download a single FASTA from NCBI Entrez."""
    url = f"{EFETCH}?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    for attempt in range(retries):
        data = _curl_fetch(url)
        if data and data.startswith(">"):
            return data
        if attempt < retries - 1:
            time.sleep(2 * (attempt + 1))
    return None


def search_ncbi(query, max_results=100):
    """Search NCBI Nucleotide and return accession list."""
    import urllib.parse
    params = urllib.parse.urlencode({
        "db": "nuccore",
        "term": query,
        "retmax": max_results,
        "retmode": "xml",
        "usehistory": "n",
    })
    url = f"{ESEARCH}?{params}"
    xml_data = _curl_fetch(url, timeout=30)
    if xml_data:
        try:
            root = ET.fromstring(xml_data)
            id_list = root.find("IdList")
            if id_list is not None:
                return [id_elem.text for id_elem in id_list.findall("Id")]
        except ET.ParseError as e:
            print(f"  XML parse failed: {e}")
    return []


def fetch_fasta_batch(gi_list, batch_size=10):
    """Fetch FASTA for a list of GI/accession IDs in batches."""
    all_fastas = {}
    for i in range(0, len(gi_list), batch_size):
        batch = gi_list[i:i + batch_size]
        ids = ",".join(batch)
        url = f"{EFETCH}?db=nuccore&id={ids}&rettype=fasta&retmode=text"
        data = _curl_fetch(url, timeout=60)
        if data:
            current_header = None
            current_seq = []
            for line in data.split("\n"):
                if line.startswith(">"):
                    if current_header:
                        acc = current_header.split()[0].lstrip(">")
                        all_fastas[acc] = f"{current_header}\n{''.join(current_seq)}\n"
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_header:
                acc = current_header.split()[0].lstrip(">")
                all_fastas[acc] = f"{current_header}\n{''.join(current_seq)}\n"
        time.sleep(0.5)  # NCBI rate limit
    return all_fastas


def main():
    print("=" * 70)
    print("Downloading Gram-positive plasmid training sequences")
    print("=" * 70)

    total_downloaded = 0

    for group_name, info in GRAM_POSITIVE_GROUPS.items():
        group_dir = os.path.join(TRAINING_DIR, group_name, "fastas")
        os.makedirs(group_dir, exist_ok=True)

        existing = [f for f in os.listdir(group_dir) if f.endswith((".fasta", ".fa"))]
        print(f"\n{'─' * 60}")
        print(f"  {group_name}: {info['description']}")
        print(f"  Organism: {info['organism']}")
        print(f"  Existing sequences: {len(existing)}")

        if len(existing) >= info["max_download"]:
            print(f"  Already have {len(existing)} sequences, skipping.")
            total_downloaded += len(existing)
            continue

        # Step 1: Download known reference accessions
        downloaded = set(f.replace(".fasta", "") for f in existing)
        for acc in info["known_accessions"]:
            if acc in downloaded:
                continue
            print(f"  Fetching known reference: {acc} ...", end=" ")
            fasta = fetch_fasta(acc)
            if fasta and len(fasta) > 100:
                out_path = os.path.join(group_dir, f"{acc}.fasta")
                with open(out_path, "w") as f:
                    f.write(fasta)
                downloaded.add(acc)
                print("OK")
            else:
                print("FAILED")
            time.sleep(0.4)

        # Step 2: Search NCBI for additional sequences
        remaining = info["max_download"] - len(downloaded)
        if remaining > 0:
            print(f"  Searching NCBI for up to {remaining} more sequences ...")
            gi_ids = search_ncbi(info["search_query"], max_results=remaining + 20)
            if gi_ids:
                print(f"  Found {len(gi_ids)} candidates, fetching ...")
                fastas = fetch_fasta_batch(gi_ids[:remaining + 10])
                for acc, fasta_text in fastas.items():
                    if acc in downloaded:
                        continue
                    if len(downloaded) >= info["max_download"]:
                        break
                    # Basic validation: must be >500 bp
                    seq_lines = [l for l in fasta_text.split("\n")
                                 if not l.startswith(">") and l.strip()]
                    seq_len = sum(len(l.strip()) for l in seq_lines)
                    if seq_len < 500:
                        continue
                    out_path = os.path.join(group_dir, f"{acc}.fasta")
                    with open(out_path, "w") as f:
                        f.write(fasta_text)
                    downloaded.add(acc)

        final_count = len([f for f in os.listdir(group_dir) if f.endswith((".fasta", ".fa"))])
        print(f"  Final count: {final_count} sequences")
        total_downloaded += final_count

    print(f"\n{'=' * 70}")
    print(f"Total Gram-positive training sequences: {total_downloaded}")
    print(f"Training directory: {TRAINING_DIR}")
    print(f"\nNext steps:")
    print(f"  1. python build_inc_centroids.py  (retrain classifier)")
    print(f"  2. streamlit run plin_app.py       (verify expanded groups)")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
