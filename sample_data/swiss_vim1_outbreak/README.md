# Sample data: Swiss VIM-1 outbreak cluster

Eight clinical isolates from a Swiss investigation of VIM-1-producing *Enterobacter* (Seth-Smith et al., Antimicrob. Agents Chemother. 70:e01827-25, 2026; reads in ENA project PRJEB98563), provided so you can try pLIN without your own FASTA files. The plasmid contigs come from Flye assemblies of the published Nanopore reads.

## What's here

- `NARACHVIM11_plasmids.fasta` ... `NARACHVIM56_plasmids.fasta`: 8 isolates, each with that isolate's assembled plasmid contigs only (19 contigs in total; chromosomal sequence removed).
- `expected_pLIN_results.tsv`: the pLIN v4.1 codes (database db-2026.10.05) produced when all 8 files are analysed together, the most similar database plasmid for each contig, and, for the blaVIM-1 plasmids, the lineage found by whole-plasmid alignment.

## How to use

1. In the app, upload all 8 `*_plasmids.fasta` files at once, keep the code scheme on **v4.1** and click **Run Analysis**.
2. Compare the **Results** tab with `expected_pLIN_results.tsv`. Codes with a value in `provisional_from` contain levels that are new to the database; their numbers are the same as in the file only when the same 8 files are analysed together in one run.

## What it shows

Every isolate carries blaVIM-1 on a large plasmid (249 to 342 kb). Whole-plasmid alignment separates them into three lineages:

| Lineage (alignment) | Isolates | pLIN v4.1 |
|---|---|---|
| A (blaVIM-1, blaCTX-M-9, blaSHV-12) | 11, 12, 20, 31, 52 | share L1 to L5 (`169.178.183.208.209`) |
| B (blaVIM-1, qnrB2) | 36, 48 | share the backbone L1 to L4 with A, own lineage at L5 |
| C (blaVIM-1, qnrB4) | 56 | same backbone family (L1 `169`), separate below |

pLIN agrees with the alignment on 27 of 28 pairs of blaVIM-1 plasmids. The one difference is isolates 11 and 52, which share 79% of their sequence at 99.997% identity, just under the 80% used to define a lineage; pLIN places them in the same lineage. Near-identical pairs (12 and 52; 20 and 31; 36 and 48) also share L6.

This is a development case study: it was used while designing pLIN v4.1, so it illustrates the method rather than testing it independently.
