# Pre-registration: protein-family coarse levels for pLIN (scheme v4)

Written 2026-10-02, before any v4 code was built or evaluated. Nothing below
may be changed after evaluation starts; any deviation is listed in a
"Deviations" section appended at the end, with the date and reason.

## 1. Question

Do pLIN codes whose coarse levels (L1–L4) are built from shared protein
families, with L5–L6 kept on 4-mer composition, group related plasmids at
least as accurately as MOB-suite and pling, while keeping pLIN's label
stability and speed?

## 2. Evidence that motivates it (pilot, already run)

`pilot_protein_backbone.py`, 2,747 plasmids, 3,183 alignment-checked pairs
(`output/protein_pilot/summary.json`):

| related-backbone AUC | all pairs | pairs split at L1–L4 |
|---|---|---|
| 4-mer cosine distance | 0.937 | 0.817 |
| backbone protein-family containment | 0.980 | 0.943 |
| all protein-family containment | 0.982 | — |

Known pilot weaknesses, which this design removes: thresholds were tuned and
scored on the same pairs, and the pairs were sampled by pLIN code sharing.

## 3. v4 scheme (fixed before evaluation)

- Proteins: Prodigal (`-p meta`), all plasmids in the release database.
- Protein families: MMseqs2 clustering at >= 50% identity, >= 80% coverage
  (pilot settings). The family catalogue is frozen per release; new plasmids'
  proteins are assigned to existing families by MMseqs2 search at the same
  thresholds, otherwise they found new families. Family IDs are never reused.
- Plasmid similarity for L1–L4: containment = shared families / families of
  the plasmid with fewer families. **All proteins** (primary). Cargo-masked
  backbone families are a secondary analysis on plasmids that have AMR/IS
  annotation (~5,800), because annotation does not exist database-wide and the
  pilot showed no loss from using all proteins.
- Assignment: the existing founder rule (`plin_founder.py`), unchanged in
  logic; at L1–L4 a plasmid joins the earliest founder with containment >=
  threshold, at L5–L6 the earliest founder within cosine distance threshold,
  always inside the cluster chosen at the level above.
- Thresholds: L1 = 0.20, L2 = 0.35, L4 = 0.75 (fixed now); L3 ("backbone
  group") = value from the calibration grid {0.40, 0.45, ..., 0.70} that
  maximises weighted F1 against "related backbone" on the calibration half
  only. L5 = 0.010, L6 = 0.001 cosine (unchanged from v3).
- v4 is a new nomenclature version. v3 codes are not silently changed; a
  v3 -> v4 mapping table is published.

## 4. Evaluation data (independent of every tool)

- 2,000 plasmids drawn uniformly at random (seed 2026) from the release
  database, excluding the 2,747 earlier benchmark plasmids. Split at random by
  **plasmid** into a calibration half and a test half (1,000 each); no plasmid
  is in both.
- Pairs: within each half, all pairs are scored with sourmash MinHash Jaccard
  (k = 21, scaled = 1000) — a measure none of the compared tools uses for
  final labels. Pairs are stratified into Jaccard bins [0, 0.01), [0.01, 0.05),
  [0.05, 0.15), [0.15, 0.4), [0.4, 0.8), [0.8, 1]; up to 600 pairs per bin per
  half are sampled at random and aligned with the existing pipeline
  (`validate_alignment_backbone.py` metrics).
- Truth (unchanged from v3 work): same lineage = bidirectional aligned
  fraction >= 0.8 and identity >= 99%; related backbone = >= 50% of the
  shorter plasmid's backbone aligns (aligned fraction of the shorter plasmid
  where backbone masking is unavailable).
- Secondary external truth: published plasmid taxonomic units (PTUs,
  Redondo-Salvo et al. 2020) for plasmids that have one, if the assignments
  can be obtained; adjusted Rand index and pairwise F1.

## 5. Methods compared (all run on the same 2,000 plasmids)

pLIN v3 (L3–L6), pLIN v4 (L1–L6), MOB-suite mob_typer primary and secondary
clusters, pling communities and subcommunities (`--sourmash`, default
thresholds). Combined rule: pLIN L6 AND MOB-suite secondary.

## 6. Endpoints

Primary (test half only; pairwise metrics weighted by inverse sampling
probability per Jaccard bin; 95% CIs from 1,000 bootstrap resamples of
plasmids):

1. Related backbone: weighted F1 of v4 L3 vs MOB-suite primary cluster.
2. Same lineage: weighted F1 of v4 L6 vs MOB-suite secondary cluster.

Secondary: unweighted metrics; precision and recall separately; per-level
results; PTU agreement; label stability (design of `benchmark_stability.py`,
5 seeds); speed (design of `benchmark_speed.py`, same machine, same threads);
share of plasmids assigned.

## 7. Decision rules

- **Gate 1 (proceed with v4):** for endpoint 1, the lower 95% bound of
  (F1 v4 L3 − F1 MOB-suite primary) is > −0.02 (non-inferior), AND v4
  incremental stability is 100% at all levels.
- **Gate 2 (Nature Communications-level claim):** Gate 1 passes, AND v4 is
  non-inferior on endpoint 2, AND v4 is faster per plasmid than MOB-suite and
  pling.
- If Gate 1 fails: v4 is reported as a negative result in the supplement, and
  the paper uses the complementary framing (pLIN as a fast, permanent naming
  layer alongside MOB-suite).

## 8. Reporting

Every number in the manuscript is produced by a script into a results file;
a checker script compares the manuscript text against those files before
submission. Results that disfavour pLIN are reported with the same
prominence as those that favour it.

## Deviations

1. 2026-10-02, before any evaluation. Case not covered above: plasmids with
   no predicted protein cannot be compared by containment. They receive
   L1–L4 = 0 and are coded at L5–L6 by 4-mer composition within that bucket.
   Their number is reported. Reason: implementation detail; no data on any
   tool's performance had been seen.
2. 2026-10-02, before any evaluation. Proteins are predicted with pyrodigal
   3.7.0 (metagenomic mode) for every plasmid instead of the Prodigal binary.
   Reason: the binary crashed (segmentation fault) on 2 of 30 randomly chosen
   database plasmids, first seen on AF318175.1. On the 28 it handled,
   pyrodigal gave identical proteins for 21 and exactly one extra protein for
   7. A single predictor for all plasmids avoids mixing two, and pyrodigal is
   what the app uses for query plasmids. No data on any tool's performance
   had been seen.
3. 2026-10-02, AFTER the first evaluation had been run and seen. While
   implementing query mode, re-querying 200 database plasmids from sequence
   reproduced only 87.5% of their L6 codes, because MMseqs2 had split
   10,309 identical protein sequences (218,604 protein copies, 1.8%) across
   different families. Every identical sequence is now placed in its lowest
   family ID (`plin_v4.py canonicalize`), the codes were rebuilt and the
   evaluation was re-run with the unchanged rules; re-querying then
   reproduces 100% of codes. Both evaluations are kept and reported
   (`evaluation_raw/` and `evaluation/`): L3 = 0.40 in both; endpoint 1
   difference +0.142 [+0.035, +0.346] before and +0.152 [+0.038, +0.352]
   after; endpoint 2 not non-inferior in either. Reason: correctness of the
   nomenclature (identical proteins must share a family), not performance.
4. 2026-10-02, after the first evaluation but before mge-cluster was run.
   mge-cluster (Arredondo-Alonso et al. 2023) is added as a secondary
   comparator: accuracy on the test half (model built on all 2,000
   evaluation plasmids, default settings) and label stability (model on
   snapshot A, then `--existing` for all plasmids, and a rebuilt model).
   Reason: a literature check found it is the closest published method
   offering consistent labels for new plasmids. It does not enter the
   primary endpoints or the gates.
