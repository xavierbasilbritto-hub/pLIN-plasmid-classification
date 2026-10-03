# Pre-registration: pLIN scheme v4.1 — confirmatory evaluation

Registered 2026-10-03, reviewed by the authors, before the confirmatory data
(§5) were drawn; no v4.1 result on those data exists. Nothing below may
change; any deviation is appended under "Deviations" with its date and reason.

## 1. Why v4.1

The pre-registered v4 evaluation (`output/backbone_v4/PREREGISTRATION.md`,
deviations 1–4) passed Gate 1, but later checks found design flaws:

- **Swiss VIM-1 outbreak:** near-identical blaVIM-1 plasmids (99.997% identity)
  were split at L1, because L1 joined the earliest founder sharing ≥ 20% of
  protein families, so weak links through widely shared mobile-element proteins
  decided the top level. Sub-lineages that alignment separates were merged at
  L6, because 4-mer composition cannot see 20–30% differences in gene content.
- **Point mutations:** after 0.1% random substitutions, only 82% of plasmids
  kept their L4 code and 80% kept L6, even with nearest-catalogued-protein
  family assignment.
- **Lineage recall:** L6 found 44% of same-lineage pairs.

v4 is therefore not released. v4.1 was designed on development data only (§3)
and is tested here on data never used before (§5).

## 2. v4.1 scheme (fixed)

Inputs per plasmid: protein families (pyrodigal 3.7.0 metagenomic mode;
family = catalogue family of the identical sequence, else of the nearest
catalogued protein by MMseqs2 search at ≥ 50% identity, ≥ 80% coverage, else
a new family; catalogue = db-2026.10.03 with identical sequences in one
family) and a k-mer sketch (`plin_kmers.py`: canonical 21-mers, splitmix64
hash, FracMinHash with a power-of-two scale chosen so ≥ 400 k-mers are kept,
maximum scale 256; two sketches are compared at the coarser scale).

Similarities
- `prot`  shared protein families / families of the plasmid with fewer
  families, counted only if ≥ min(3, that count) families are shared
- `kcont` k-mer containment of the plasmid with fewer k-mers in the other
- `kmin`  shared k-mers / k-mers of the larger plasmid (symmetric)

Levels and assignment (`plin_v41_nn.HybridIndex`)
| level | meaning | similarity | threshold | rule |
|---|---|---|---|---|
| L1 | backbone family | prot | 0.40 | founder (earliest) |
| L2 | backbone group | prot | 0.60 | founder (earliest) |
| L3 | shared backbone | kcont | 0.50 | founder (earliest) |
| L4 | backbone variant | kcont | 0.80 | founder (earliest) |
| L5 | lineage | kmin | 0.80 | nearest-neighbour copy |
| L6 | near-identical (outbreak clone) | kmin | 0.95 | nearest-neighbour copy |

A plasmid whose most similar database plasmid (kmin; ties → earliest) reaches
0.80 copies that plasmid's L1–L5 code, then copies L6 from the most similar
plasmid inside that lineage if kmin ≥ 0.95, else gets a new L6 ID. Otherwise
it is placed at L1–L4 by the founder rule (each level searches only inside
the cluster chosen above; plasmids without proteins get L1 = L2 = 0) and gets
new L5 and L6 IDs. The release database is built in accession order;
database plasmids keep their codes in later releases.

## 3. Development (already done; not evidence of performance)

Data used during design, all excluded from §5: the 2,747 comparator-benchmark
plasmids, the 2,000 v4 evaluation plasmids and their 5,544 aligned pairs,
50 mutation-test plasmids, 200 re-query-check plasmids, the 74 published
outbreak plasmids and the 19 Swiss VIM-1 plasmid contigs. Designs compared
(`output/backbone_v41/dev/compare.tsv`): founder-only (best / earliest
founder), nearest-neighbour-only (6 threshold sets) and hybrid (6 threshold
sets). The design in §2 (H3) was chosen for the best combination of backbone
F1, lineage F1, re-query reproduction and mutation robustness. Swiss VIM-1
and the outbreak set are reported as development case studies, not as
confirmation.

## 4. Release database

db-2026.10.03: 127,517 unique plasmids (5,788 accessions present twice in
db-2026.10.02, as X and RefSeq_X, are merged).

## 5. Confirmatory data

- 2,000 plasmids drawn uniformly at random (seed 4141) from the release,
  excluding every plasmid listed in §3 and any plasmid whose accession matches
  them.
- Truth pairs: the v4 design (`v4_truth_pairs.py`): all pairs scored with
  sourmash MinHash Jaccard (k = 21, scaled = 1000), stratified into bins
  [0, 0.01), [0.01, 0.05), [0.05, 0.15), [0.15, 0.4), [0.4, 0.8), [0.8, 1],
  up to 600 pairs per bin, aligned with blastn; weights = inverse sampling
  probability. Same lineage: AF_min ≥ 0.8 and identity ≥ 99%. Related
  backbone: ≥ 50% of the shorter plasmid aligned (backbone-masked where
  annotated). All 2,000 plasmids form one test set; nothing is tuned on it.
- Mutation set: 100 further fresh plasmids (seed 4242), random substitutions
  at 0.01%, 0.1% and 1% (seed 4343).

## 6. Comparators (same 2,000 plasmids)

pLIN v4 (pre-registered baseline), MOB-suite mob_typer (primary and secondary
clusters), pling (community, subcommunity; `--sourmash`, defaults),
mge-cluster 1.1.0 (defaults, model on the 2,000).

## 7. Endpoints

Primary (weighted; 95% CIs from 1,000 bootstrap resamples of plasmids):
1. Related backbone: F1 of v4.1 L1 vs MOB-suite primary cluster.
2. Same lineage: F1 of v4.1 L5 vs MOB-suite secondary cluster.

Secondary: unweighted metrics; precision and recall; all levels and all
comparators; agreement with published PTUs where available; outbreak
false-link rate per level is reported only descriptively (development set).

## 8. Release criteria (all must hold; checked on the confirmatory data)

- R1 Stability: 100% of codes unchanged at every level when the database grows
  (snapshots 25/50/75%, seeds 1–5).
- R2 Reproduction: 100% of 500 fresh release plasmids re-query to their own
  code through the full sequence path (no whole-plasmid shortcut).
- R3 Robustness: ≥ 95% of mutated plasmids keep their original L1–L5 code at
  0.01% and at 0.1% substitutions (1% reported only).
- R4 Integrity: no duplicate accessions; every code reproduced by the compact
  release files; SHA-256 checksums published.
- R5 Accuracy: primary endpoint 1 non-inferior to MOB-suite primary (lower
  95% bound of the F1 difference > −0.02).

If R1–R4 fail, v4.1 is not released. If R5 fails, or endpoint 2 is not
non-inferior to MOB-suite secondary, the result is reported as a limitation
and the paper does not claim superiority on that endpoint.

## Deviations

(none yet)
