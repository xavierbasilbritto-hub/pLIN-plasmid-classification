#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0-or-later
# See LICENSE. Please cite pLIN (CITATION.cff).
"""
One Health view of pLIN v4 lineages (descriptive; not pre-registered).

Release plasmids are matched to PLSDB 2024_05_31_v2 metadata (figshare
10.6084/m9.figshare.27252609) and given a One Health sector:
  human        BioSample host taxid 9606
  animal       host taxid within Metazoa (33208), other than human
  plant        host taxid within Viridiplantae (33090)
  food         no host; ecosystem tags include food, meat, fermented or drink
  environment  no host; tags include soil, aquatic, freshwater, marine, sea,
               river, lake, sediment, wastewater, sewage, wwtp, terrestrial,
               saline, contaminated or hospital (built environment)
  unknown      anything else
Host lineages come from NCBI Taxonomy (E-utilities), cached in
taxid_lineage.json. AMR genes are PLSDB's AMRFinderPlus calls, grouped as
carbapenemase, mcr (colistin), CTX-M (ESBL) and tet(X) (tigecycline).

For v4 L6 lineages (and L3 backbone groups) with >= 2 annotated plasmids:
how many span >= 2 sectors, and how many carry the same key AMR gene in >= 2
sectors. A shared code across sectors is consistent with spread of the same
plasmid lineage but does not by itself establish transmission (no dates,
sampling bias towards clinical isolates).

Usage:
  python v4_onehealth.py
Output: output/backbone_v4/onehealth/{plasmid_sectors.tsv, lineage_summary.tsv,
        cross_sector_examples.tsv, summary.json}
"""

import json
import os
import re
import time
import urllib.request
import xml.etree.ElementTree as ET

import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OH = os.environ.get("PLIN_ONEHEALTH_OUT", os.path.join(BASE_DIR, "output", "backbone_v4", "onehealth"))
META = os.path.join(BASE_DIR, "output", "backbone_v4", "onehealth", "plsdb_meta")
CODES = os.environ.get("PLIN_ONEHEALTH_CODES",
                       os.path.join(BASE_DIR, "output", "backbone_v4", "codes", "plin_v4_codes_L3_0.40.tsv"))
CODE_COL = os.environ.get("PLIN_ONEHEALTH_CODE_COL", "pLIN_v4")
ENV_TAGS = {"soil", "aquatic", "freshwater", "marine", "sea", "river", "lake", "sediment", "wastewater",
            "sewage", "wwtp", "terrestrial", "saline", "contaminated", "hospital"}
FOOD_TAGS = {"food", "meat", "fermented", "drink"}
AMR_GROUPS = {"carbapenemase": r"^bla(?:KPC|NDM|VIM|IMP|GES-(?:5|14|16|18|20)|OXA-(?:48|162|163|181|204|232|244|245|484|515|535))\b",
              "mcr": r"^mcr-", "CTX-M": r"^blaCTX-M", "tet(X)": r"^tet\(X"}


def taxid_lineages(taxids):
    """taxid -> list of ancestor taxids (NCBI E-utilities), cached."""
    cache_f = os.path.join(OH, "taxid_lineage.json")
    cache = json.load(open(cache_f)) if os.path.exists(cache_f) else {}
    todo = [t for t in taxids if str(t) not in cache]
    for k in range(0, len(todo), 150):
        ids = ",".join(map(str, todo[k:k + 150]))
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={ids}&retmode=xml"
        root = ET.fromstring(urllib.request.urlopen(url, timeout=60).read())
        for tx in root.findall("Taxon"):
            tid = tx.findtext("TaxId")
            cache[tid] = [x.findtext("TaxId") for x in tx.findall("LineageEx/Taxon")] + [tid]
        time.sleep(0.4)                                    # NCBI rate limit without an API key
    json.dump(cache, open(cache_f, "w"))
    return cache


def sectors():
    n = pd.read_csv(os.path.join(META, "nuccore.csv"), usecols=["NUCCORE_ACC", "BIOSAMPLE_UID"])
    b = pd.read_csv(os.path.join(META, "biosample.csv"), low_memory=False,
                    usecols=["BIOSAMPLE_UID", "ECOSYSTEM_taxid", "ECOSYSTEM_tags", "LOCATION_name"]) \
        .drop_duplicates("BIOSAMPLE_UID")
    d = n.merge(b, on="BIOSAMPLE_UID", how="left")
    d["taxid"] = pd.to_numeric(d.ECOSYSTEM_taxid, errors="coerce").fillna(-1).astype(int)
    lin = taxid_lineages(sorted(t for t in d.taxid.unique() if t > 0))

    def sector(taxid, tags):
        tags = set(str(tags).split(",")) if pd.notna(tags) else set()
        if taxid == 9606:
            return "human"
        if taxid > 0:
            anc = lin.get(str(taxid), [])
            return "animal" if "33208" in anc else "plant" if "33090" in anc else "unknown"
        if tags & FOOD_TAGS:
            return "food"
        if tags & ENV_TAGS:
            return "environment"
        return "unknown"
    d["sector"] = [sector(t, g) for t, g in zip(d.taxid, d.ECOSYSTEM_tags)]
    return d[["NUCCORE_ACC", "sector", "LOCATION_name"]]


def amr_groups():
    a = pd.read_csv(os.path.join(META, "amr.tsv"), sep="\t", usecols=["NUCCORE_ACC", "gene_symbol"])
    rows = []
    for g, rx in AMR_GROUPS.items():
        hit = a[a.gene_symbol.astype(str).str.contains(rx, regex=True)]
        rows.append(hit.assign(group=g)[["NUCCORE_ACC", "group", "gene_symbol"]])
    return pd.concat(rows).drop_duplicates()


def main():
    codes = pd.read_csv(CODES, sep="\t").rename(columns={CODE_COL: "pLIN_v4"})
    codes["acc"] = codes.plasmid_id.str.replace("^RefSeq_", "", regex=True)
    # 5,788 accessions are in the release twice (with and without "RefSeq_"; same sequence, same code)
    codes = codes.sort_values("plasmid_id").drop_duplicates("acc")
    sec = sectors()
    d = codes.merge(sec, left_on="acc", right_on="NUCCORE_ACC", how="inner")
    amr = amr_groups()
    d.to_csv(os.path.join(OH, "plasmid_sectors.tsv"), sep="\t", index=False)
    known = d[d.sector != "unknown"]
    summary = {"release_unique_accessions": len(codes), "matched_to_PLSDB_2024_05_31": len(d),
               "sector_counts": d.sector.value_counts().to_dict()}

    examples = []
    for level in (6, 3):
        known = known.assign(lineage=known.pLIN_v4.map(lambda c, k=level: ".".join(c.split(".")[:k])))
        # drop the no-protein bucket at protein levels
        if level <= 4:
            known = known[~known.lineage.str.startswith("0")]
        g = known.groupby("lineage")
        multi = g.filter(lambda x: len(x) >= 2)
        span = multi.groupby("lineage").sector.nunique()
        ka = multi.merge(amr, left_on="acc", right_on="NUCCORE_ACC")
        gene_span = ka.groupby(["lineage", "gene_symbol"]).sector.nunique()
        cross_gene = gene_span[gene_span >= 2].reset_index()
        summary[f"L{level}"] = {"lineages_with_>=2_annotated_plasmids": int(span.size),
                                "spanning_>=2_sectors": int((span >= 2).sum()),
                                "same_key_AMR_gene_in_>=2_sectors": int(cross_gene.lineage.nunique())}
        for (lineage, gene), _ in cross_gene.set_index(["lineage", "gene_symbol"]).iterrows():
            m = ka[(ka.lineage == lineage) & (ka.gene_symbol == gene)]
            examples.append({"level": f"L{level}", "lineage": lineage, "gene": gene,
                             "group": m.group.iloc[0], "plasmids_with_gene": m.acc.nunique(),
                             "sectors": ",".join(sorted(m.sector.unique())),
                             "countries": m.LOCATION_name.dropna().nunique(),
                             "lineage_size": int(g.size()[lineage])})
    ex = pd.DataFrame(examples).sort_values(["level", "plasmids_with_gene"], ascending=[False, False])
    ex.to_csv(os.path.join(OH, "cross_sector_examples.tsv"), sep="\t", index=False)
    json.dump(summary, open(os.path.join(OH, "summary.json"), "w"), indent=2)
    print(json.dumps(summary, indent=2))
    pd.set_option("display.width", 200)
    print(ex.groupby("level").head(10).to_string(index=False))


if __name__ == "__main__":
    main()
