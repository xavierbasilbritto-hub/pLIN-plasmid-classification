# Table 1. Comparative evaluation of plasmid classification systems against five criteria required for comprehensive plasmid nomenclature.

| Property | pLIN (this study) | PlasmidFinder / Inc typing [9] | pMLST [9] | MOB-suite [10] | COPLA / PTUs [11] | mge-cluster [12] |
|---|---|---|---|---|---|---|
| **Year introduced** | 2025 | 2014 | 2014 | 2018 | 2021 | 2023 |
| **Classification basis** | Whole-sequence 4-mer composition (256 features) | Replicon gene detection (BLAST) | Allelic variants of replicon loci | Relaxase typing + Mash distance | ANI network + hierarchical stochastic block modelling | Unitig Jaccard distance + HDBSCAN |
| **Hierarchical levels** | 6 nested levels (L1--L6) | 1 (flat) | 1 (flat) | 1 (flat) | 2--3 (semi-hierarchical) | 1 (flat) |
| **Code permanence** | **Yes** (guaranteed by nearest-neighbour rule) | Stable (DB-dependent) | Stable (DB-dependent) | **No** (re-clustering on DB update) | **No** (HSBM recomputation) | **No** (t-SNE re-embedding) |
| **Reference database required** | **No** (reference-free core) | Yes | Yes | Yes | Yes | **No** |
| **Taxonomic scope** | 28 groups validated (20 GN Inc + 4 GP rep + 2 Aci + 2 Pae); universal (any plasmid) | Enterobacteriaceae mainly (~30 replicons) | 6 Inc schemes only | Broad (DB-limited) | Broad (41% assignable; 63% Enterobacterales) | Broad |
| **Resolution** | Strain-level (~99.9% ANI) | Family-level | Sub-family (ST) | Species-level (Mash 0.06) | Species-level | Variable |
| **Simpson's D (this dataset, n = 8,077)** | **0.985** | 0.641 | N/A^a^ | N/A | N/A | N/A |
| **Improvement over Inc/rep typing** | 1.54-fold | -- | -- | -- | -- | -- |
| **Sensitivity** | High (whole-sequence) | ~50% (replicon-dependent) | Scheme-limited | ~95% | 41--63% | High |
| **Computational cost** | Low (minutes) | Low (seconds) | Low (seconds) | Moderate (minutes) | High (hours) | Moderate (minutes) |
| **ML validation** | Integrated (XGBoost F1 = 0.896) | None | None | None | None | None |
| **AMR integration** | **Yes** (AMRFinderPlus) | No | No | No | No | No |
| **Outbreak detection** | **Yes** (basic + temporal) | No | No | No | No | No |
| **Mosaic detection** | **Yes** (GC heterogeneity) | No | No | No | No | No |

^a^ pMLST schemes not available for all 20 Gram-negative Inc groups simultaneously.

**Notes:**

- Simpson's Index of Diversity (D) calculated on the full 8,077-sample, 28-group training dataset. Inc/rep group typing produces D = 0.641 because it assigns all plasmids within each group to a single label; pLIN resolves 3,073 unique strain-level codes across all 28 groups.
- Code permanence for PlasmidFinder/pMLST is labelled "Stable (DB-dependent)" because the codes themselves do not change, but classification outcomes can change when the underlying replicon/allele database is updated.
- The "reference-free" designation for pLIN refers to the core classification algorithm, which operates solely on k-mer composition without requiring external databases. The optional AMRFinderPlus, MOB-suite, and FastANI modules do require their respective databases.

**Abbreviations:** ANI, average nucleotide identity; DB, database; HDBSCAN, hierarchical density-based spatial clustering of applications with noise; HSBM, hierarchical stochastic block modelling; ML, machine learning; pMLST, plasmid multilocus sequence typing; PTU, plasmid taxonomic unit; ST, sequence type.
