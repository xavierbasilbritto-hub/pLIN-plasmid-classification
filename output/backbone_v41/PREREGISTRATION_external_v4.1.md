# Pre-registration: pLIN v4.1, external validation on independent hospital datasets

Registered 2026-10-04, before any sequence of the datasets below was downloaded or typed. The method
(pLIN v4.1, database db-2026.10.03) and its thresholds are frozen as registered in
PREREGISTRATION_v4.1.md. Nothing below may change; deviations are appended with date and reason.

## 1. Purpose

The outbreak analyses of the confirmatory study used data that informed the design. This analysis tests
pLIN v4.1 on published hospital datasets never used during development.

## 2. Datasets

- **E1**: multispecies hospital outbreak of NDM-5-producing Enterobacterales, USA, 2021 to 2023
  (hybrid assemblies, NCBI BioProject PRJNA981541).
- **E2**: endemic IMP-4 carbapenemase dissemination through clonal, plasmid and integron transfer in a
  hospital network, Melbourne, Australia (Macesic et al., Nat Commun 2023; NCBI BioProject PRJNA924056).

Inclusion: all nuccore records linked to the BioProject (directly or through its assemblies) whose title
contains "plasmid", that are complete sequences, and that are at least 1,000 bp. Excluded: any accession
(with or without version) in the development data of PREREGISTRATION_v4.1.md section 3 or in the
confirmatory test, mutation or re-query sets. If a BioProject yields fewer than 20 plasmids, it is
reported but not analysed.

## 3. Typing

- **pLIN v4.1**: all plasmids of a dataset typed in one session from sequence (plin_v41_typer), with every
  database plasmid of identical sequence (whole-sequence hash) or identical accession treated as absent,
  as for a hospital receiving new isolates.
- **MOB-suite** mob_typer 3.1.9, default settings (primary and secondary clusters).
- **pling** 3.0.2, default settings with --sourmash, run per dataset (community and subcommunity).

## 4. Truth

As registered for the confirmatory study (v4_truth_pairs.py design), within each dataset: all pairs scored
by sourmash MinHash Jaccard (k = 21, scaled = 1,000), stratified into the same six Jaccard bins, up to 600
pairs per bin aligned with blastn; weights = inverse sampling probability.
Same lineage: AF_min >= 0.8 and identity >= 99%. Related backbone: >= 50% of the shorter plasmid
aligned (backbone-masked where annotated).

## 5. Endpoints

Primary (pooled over E1 and E2; 95% CIs from 1,000 bootstrap resamples of plasmids within datasets):
1. Same lineage: weighted F1 of pLIN L5 minus that of MOB-suite secondary clusters.
2. Related backbone: weighted F1 of pLIN L1 minus that of MOB-suite primary clusters.

Secondary: per-dataset F1, precision and recall for every method and level; for E1, the share of pairs of
NDM-5-carrying plasmids that share L5 and L6 (the authors describe one outbreak plasmid); for E2, the
same for IMP-4-carrying plasmids, by the authors' plasmid types where available.

All results are reported whatever their direction.

## Deviations

(none yet)
