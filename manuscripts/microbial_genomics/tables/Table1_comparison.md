# Table 1. Comparative evaluation of plasmid classification and typing systems.

| Property | pLIN (this study) | PlasmidFinder / Inc typing | pMLST | MOB-suite | COPLA (PTUs) | mge-cluster |
|---|---|---|---|---|---|---|
| **Year introduced** | 2025 | 2014 | 2014 | 2018 | 2020 | 2023 |
| **Classification basis** | Whole-sequence tetranucleotide (4-mer) composition + cosine distance | Replicon gene detection (BLAST) | Allelic variants of replicon loci | Relaxase typing + Mash distance | ANI network + hierarchical stochastic block modelling | Unitig Jaccard distance + HDBSCAN |
| **Hierarchical levels** | 6 nested levels (L1--L6) | 1 (flat) | 1 (flat) | 1 (flat) | 2--3 (semi-hierarchical) | 1 (flat) |
| **Code permanence** | Yes (guaranteed by nearest-neighbour rule) | Yes (database-dependent) | Yes (database-dependent) | No (re-clustering on update) | No (HSBM recomputation) | No (re-embedding required) |
| **Reference database required** | No | Yes | Yes | Yes | Yes | No |
| **Taxonomic scope** | 28 groups validated (20 GN Inc + 4 GP rep + 2 Aci + 2 Pae); universal in principle | Enterobacteriaceae mainly | 6 Inc group schemes only | Broad (database-limited) | Broad (41% assignable) | Broad |
| **Resolution** | Strain-level (~99.9% ANI) | Family-level | Sub-family | Species-level (Mash d = 0.06) | Species-level | Variable |
| **Simpson's D (this dataset)** | 0.985 | 0.641* | N/A** | N/A | N/A | N/A |
| **Classification rate** | 97.3% | ~50% (replicon-dependent) | Scheme-limited | ~95% | 41--63% | High |
| **Computational cost** | Low (< 30 min for 79,305 plasmids) | Low | Low | Moderate | High (hours) | Moderate |
| **ML validation** | Integrated (XGBoost F1 = 0.896) | None | None | None | None | None |
| **AMR gene integration** | Yes (AMRFinderPlus) | No | No | No | No | No |
| **Outbreak detection** | Yes (pLIN L6 + AMR fingerprint + temporal) | No | No | No | No | No |

\* Calculated for 28-group Inc/rep typing on this dataset (n = 8,077).

\*\* pMLST schemes not available for all 20 Gram-negative Inc groups simultaneously.

PlasmidFinder and pMLST [13]; MOB-suite [14]; COPLA [15]; mge-cluster [16].
