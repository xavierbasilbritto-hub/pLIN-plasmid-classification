# Sample data: Swiss VIM-1 outbreak cluster

Eight clinical isolates from the Swiss carbapenemase-producing outbreak case study used as a flagship validation example in the pLIN manuscript, provided so you can try pLIN immediately without needing your own FASTA files.

## What's here

- `NARACHVIM11_plasmids.fasta` ... `NARACHVIM56_plasmids.fasta` — 8 isolates, each containing that isolate's assembled plasmid contigs only (chromosomal sequence already removed).
- `expected_pLIN_results.tsv` — the pLIN codes this analysis is expected to produce, from the original manuscript validation run, so you can check your own run matches.

## How to use

1. In the app's **Overview** tab, upload all 8 `*_plasmids.fasta` files at once (multi-file upload is supported).
2. Leave **Incompatibility Group** on **Auto-detect** and click **Run Analysis**.
3. Check the **Results** tab against `expected_pLIN_results.tsv` — most contigs should resolve to pLIN code **1.1.2.4.7.13**, the shared outbreak lineage (IncN/IncHI2), demonstrating that pLIN correctly identifies these epidemiologically-linked isolates as one lineage despite being typed under two different replicon families.
4. Try the **Cladogram** tab to see the isolates cluster together visually, and **Epidemiology** to see the outbreak/clone-detection module flag the shared code.

## Background

This is the same dataset behind Figure 12 (transmission mode) and the outbreak-validation section of the manuscript: a real hospital cluster where a blaVIM-1-carrying plasmid spread across patients, captured here as assembled contigs from each isolate's sequencing run.
