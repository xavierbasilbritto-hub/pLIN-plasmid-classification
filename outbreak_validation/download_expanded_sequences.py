#!/usr/bin/env python3
"""Download outbreak plasmid sequences from NCBI for expanded pLIN validation.

Accessions sourced from published outbreak studies covering 7 resistance mechanisms
across 15+ countries. Each accession is a complete plasmid sequence from a
peer-reviewed outbreak investigation.
"""

import os
import sys
import time
import urllib.request
import urllib.error

OUT_DIR = os.path.join(os.path.dirname(__file__), "expanded_sequences")
os.makedirs(OUT_DIR, exist_ok=True)

# Accessions organized by resistance mechanism and study
ACCESSIONS = {
    # KPC carbapenemases
    "GU595196": {"gene": "blaKPC-2", "country": "USA", "study": "Kitchel_2009_KPC_IncN"},
    "JN233704": {"gene": "blaKPC-2", "country": "USA", "study": "Chen_2012_KPC_IncFIA"},
    "CP004366": {"gene": "blaKPC-3", "country": "USA", "study": "Conlan_2014_NIH_KPC"},
    "CP004367": {"gene": "blaKPC-3", "country": "USA", "study": "Conlan_2014_NIH_KPC"},
    "CP019026": {"gene": "blaKPC-2", "country": "China", "study": "Li_2018_KPC_IncN"},
    "CP081510": {"gene": "blaKPC-2", "country": "Italy", "study": "Arcari_2023_KPC_outbreak"},
    "CP081509": {"gene": "blaKPC-2", "country": "Italy", "study": "Arcari_2023_KPC_outbreak"},
    "MH133192": {"gene": "blaKPC-2", "country": "Brazil", "study": "Andrade_2019_KPC_Brazil"},

    # NDM metallo-beta-lactamases (Hong Kong outbreak)
    "MH234497": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234498": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234499": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234500": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234501": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234502": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234503": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234504": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234505": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234506": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234507": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234508": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},
    "MH234509": {"gene": "blaNDM-1", "country": "Hong Kong", "study": "Ho_2019_NDM_HK_ICU"},

    # NDM - Colombia
    "KX832926": {"gene": "blaNDM-1", "country": "Colombia", "study": "Rojas_2017_NDM_Colombia"},
    "KX832927": {"gene": "blaNDM-1", "country": "Colombia", "study": "Rojas_2017_NDM_Colombia"},
    "KX832928": {"gene": "blaNDM-1", "country": "Colombia", "study": "Rojas_2017_NDM_Colombia"},
    "KX832929": {"gene": "blaNDM-1", "country": "Colombia", "study": "Rojas_2017_NDM_Colombia"},
    "CP017672": {"gene": "blaNDM-1", "country": "Colombia", "study": "Rojas_2017_NDM_Colombia"},

    # NDM - China
    "JX104760": {"gene": "blaNDM-1", "country": "China", "study": "Ho_2012_NDM_China"},
    "MH985166": {"gene": "blaNDM-5", "country": "China", "study": "Li_2020_NDM5_China"},
    "MH985167": {"gene": "blaNDM-5", "country": "China", "study": "Li_2020_NDM5_China"},
    "MH985168": {"gene": "blaNDM-5", "country": "China", "study": "Li_2020_NDM5_China"},
    "MH985169": {"gene": "blaNDM-5", "country": "China", "study": "Li_2020_NDM5_China"},
    "MH985170": {"gene": "blaNDM-5", "country": "China", "study": "Li_2020_NDM5_China"},
    "MH985171": {"gene": "blaNDM-5", "country": "China", "study": "Li_2020_NDM5_China"},

    # OXA-48
    "JN626286": {"gene": "blaOXA-48", "country": "Turkey", "study": "Potron_2013_OXA48_Turkey"},
    "LR025096": {"gene": "blaOXA-48", "country": "Netherlands", "study": "Jousset_2019_OXA48_NL"},
    "LR025097": {"gene": "blaOXA-48", "country": "Netherlands", "study": "Jousset_2019_OXA48_NL"},
    "LR025098": {"gene": "blaOXA-48", "country": "Netherlands", "study": "Jousset_2019_OXA48_NL"},
    "LR025099": {"gene": "blaOXA-48", "country": "Netherlands", "study": "Jousset_2019_OXA48_NL"},
    "LR025100": {"gene": "blaOXA-48", "country": "Netherlands", "study": "Jousset_2019_OXA48_NL"},
    "LR025105": {"gene": "blaOXA-48", "country": "Netherlands", "study": "Jousset_2019_OXA48_NL"},
    "KP061858": {"gene": "blaOXA-48", "country": "France", "study": "Jousset_2019_OXA48_FR"},

    # VIM
    "MN783743": {"gene": "blaVIM-1", "country": "Italy", "study": "Arcari_2020_VIM_Italy"},
    "MN783744": {"gene": "blaVIM-1", "country": "Italy", "study": "Arcari_2020_VIM_Italy"},
    "MN783745": {"gene": "blaVIM-1", "country": "Italy", "study": "Arcari_2020_VIM_Italy"},

    # IMP
    "AB616660": {"gene": "blaIMP-6", "country": "Japan", "study": "Tada_2015_IMP_Japan"},

    # mcr (colistin resistance)
    "KP347127": {"gene": "mcr-1", "country": "China", "study": "Liu_2016_mcr1_discovery"},
    "KU761326": {"gene": "mcr-1", "country": "China", "study": "Zheng_2017_mcr1_China"},
    "KU761327": {"gene": "mcr-1", "country": "China", "study": "Zheng_2017_mcr1_China"},
    "KY075653": {"gene": "mcr-1", "country": "Europe", "study": "Hasman_2015_mcr1_Europe"},
    "KY075654": {"gene": "mcr-1", "country": "Europe", "study": "Hasman_2015_mcr1_Europe"},
    "CP016405": {"gene": "mcr-1", "country": "USA", "study": "McGann_2016_mcr1_USA"},

    # CTX-M ESBL
    "AY458016": {"gene": "blaCTX-M-15", "country": "India", "study": "Karim_2001_CTXM_India"},
    "EU935738": {"gene": "blaCTX-M-15", "country": "UK", "study": "Woodford_2009_CTXM_UK"},
    "EU935739": {"gene": "blaCTX-M-15", "country": "UK", "study": "Woodford_2009_CTXM_UK"},
    "EU935740": {"gene": "blaCTX-M-15", "country": "UK", "study": "Woodford_2009_CTXM_UK"},
    "CP009231": {"gene": "blaCTX-M-15", "country": "USA", "study": "Sheppard_2016_CTXM_USA"},
    "CP009232": {"gene": "blaCTX-M-15", "country": "USA", "study": "Sheppard_2016_CTXM_USA"},
    "CP009233": {"gene": "blaCTX-M-15", "country": "USA", "study": "Sheppard_2016_CTXM_USA"},
    "FN868832": {"gene": "blaCTX-M-15", "country": "France", "study": "Valverde_2009_CTXM_FR"},
}


def fetch_fasta(accession: str, retries: int = 3) -> str | None:
    """Fetch a FASTA sequence from NCBI efetch."""
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nucleotide&id={accession}&rettype=fasta&retmode=text"
    )
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url)
            req.add_header("User-Agent", "pLIN_validation/1.0")
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read().decode("utf-8")
                if data.startswith(">"):
                    return data
                elif "Nothing has been found" in data or len(data.strip()) < 20:
                    return None
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            print(f"  Attempt {attempt+1} failed: {e}")
            time.sleep(2 * (attempt + 1))
    return None


def validate_sequence(fasta: str, accession: str) -> dict:
    """Basic validation of downloaded sequence."""
    lines = fasta.strip().split("\n")
    header = lines[0]
    seq = "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
    length = len(seq)

    return {
        "accession": accession,
        "header": header[:120],
        "length": length,
        "is_plasmid": length < 500000,  # < 500 kb likely plasmid
        "valid_bases": all(c in "ATCGNatcgn" for c in seq[:1000]),
    }


def main():
    total = len(ACCESSIONS)
    downloaded = 0
    failed = []
    skipped_chromo = []
    results = []

    print(f"Downloading {total} outbreak plasmid sequences from NCBI...")
    print(f"Output directory: {OUT_DIR}\n")

    for i, (acc, meta) in enumerate(ACCESSIONS.items(), 1):
        out_path = os.path.join(OUT_DIR, f"{acc}.fasta")

        # Skip if already downloaded
        if os.path.exists(out_path) and os.path.getsize(out_path) > 100:
            print(f"[{i}/{total}] {acc} — already exists, skipping")
            downloaded += 1
            continue

        print(f"[{i}/{total}] {acc} ({meta['gene']}, {meta['country']}) ... ", end="", flush=True)
        fasta = fetch_fasta(acc)

        if fasta is None:
            print("FAILED")
            failed.append(acc)
            continue

        info = validate_sequence(fasta, acc)
        results.append({**info, **meta})

        if not info["is_plasmid"]:
            print(f"SKIPPED (chromosome? {info['length']:,} bp)")
            skipped_chromo.append(acc)
            continue

        if not info["valid_bases"]:
            print(f"WARNING (invalid bases) — saving anyway")

        with open(out_path, "w") as f:
            f.write(fasta)

        downloaded += 1
        print(f"OK ({info['length']:,} bp)")

        # Rate limit: NCBI allows 3 requests/sec without API key
        time.sleep(0.4)

    print(f"\n{'='*60}")
    print(f"Downloaded: {downloaded}/{total}")
    print(f"Failed:     {len(failed)}")
    print(f"Skipped (>500kb): {len(skipped_chromo)}")

    if failed:
        print(f"\nFailed accessions: {', '.join(failed)}")
    if skipped_chromo:
        print(f"\nSkipped (likely chromosomal): {', '.join(skipped_chromo)}")

    # Write summary
    summary_path = os.path.join(OUT_DIR, "download_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("accession\tgene\tcountry\tstudy\tlength\tis_plasmid\tstatus\n")
        for r in results:
            status = "ok" if r["is_plasmid"] else "skipped_chromo"
            f.write(f"{r['accession']}\t{r.get('gene','')}\t{r.get('country','')}\t"
                    f"{r.get('study','')}\t{r['length']}\t{r['is_plasmid']}\t{status}\n")
        for acc in failed:
            meta = ACCESSIONS[acc]
            f.write(f"{acc}\t{meta['gene']}\t{meta['country']}\t{meta['study']}\t0\tFalse\tfailed\n")

    print(f"\nSummary written to: {summary_path}")
    return downloaded, failed, skipped_chromo


if __name__ == "__main__":
    main()
