# Appendix p 16: Classification coverage comparison and characterisation of unclassified sequences

## Table S6. Classification coverage comparison across plasmid typing methods

| Method | Approach | Classification coverage | Permanent codes | Hierarchical | Reference |
|--------|----------|------------------------|:-:|:-:|-----------|
| PlasmidFinder | Replicon marker detection | ~40% | Yes | No | Carattoli *et al.*, 2014 (Ref 4) |
| PlasmidFinder + MOB-typer | Replicon + relaxase detection | ~74% | Partial | No | Schmartz *et al.*, 2025 (Ref 32) |
| pMLST | Allelic profiling (6 Inc groups) | ~19% (eligible fraction) | Yes | No | Carattoli *et al.*, 2014 (Ref 4) |
| COPLA | PTU assignment (genomic taxonomy) | 41% (63% Enterobacterales) | Yes | No | Redondo-Salvo *et al.*, 2021 |
| MOB-typer | Relaxase typing | ~55% | Yes | No | Robertson & Nash, 2020 |
| mge-cluster | Reference-free clustering | 79.8% | No | No | Arredondo-Alonso *et al.*, 2023 (Ref 8) |
| MOB-cluster | Mash distance clustering | 100% | No | No | Robertson & Nash, 2020 |
| pling | Rearrangement distance | Variable (study-dependent) | No | No | Watts *et al.*, 2024 |
| **pLIN** | **4-mer cosine + hierarchical LIN** | **97.3%** | **Yes** | **Yes (6 levels)** | **This study** |

Notes: Coverage percentages for PlasmidFinder and MOB-typer are based on analysis of PLSDB (Galata *et al.*, 2019; Schmartz *et al.*, 2022). pMLST eligibility is limited to plasmids belonging to one of six typed Inc groups (IncA/C, IncHI1, IncHI2, IncI1, IncN, IncF). COPLA coverage from analysis of 10,696 plasmids (Redondo-Salvo *et al.*, 2021). mge-cluster assigns 20.2% of plasmids to noise (Arredondo-Alonso *et al.*, 2023). MOB-cluster achieves 100% assignment but cluster identifiers are regenerated with each database update, precluding longitudinal surveillance. pLIN coverage calculated as the proportion of 79,305 plasmids receiving a classification with >= 40% KNN confidence (77,196/79,305 = 97.3%).

---

## Characterisation of unclassified sequences

Of 79,305 plasmids in the expanded pLIN database, 2,109 (2.7%) received classifications below the 40% confidence threshold. These sequences were assigned forced Inc type labels and pLIN codes but are flagged as low-confidence to enable user-side quality assessment.

### Table S7. Comparison of classified versus unclassified sequences

| Property | Classified (n = 77,196) | Unclassified (n = 2,109) | Ratio |
|----------|------------------------|-------------------------|-------|
| Proportion | 97.3% | 2.7% | — |
| Median length (bp) | 59,311 | 12,516 | 4.7x shorter |
| Median NN distance | 0.0067 | 0.0206 | 3.1x more distant |
| Unique pLIN codes | 55,913 (72.4%) | 1,991 (94.4%) | — |
| Source | Training + reference | Reference only | — |

### Table S8. Nearest-neighbour distance distribution for unclassified sequences

| Statistic | Unclassified NN distance |
|-----------|-------------------------|
| Minimum | 0.0002 |
| Q1 (25th percentile) | 0.0094 |
| Median | 0.0206 |
| Q3 (75th percentile) | 0.0350 |
| Maximum | 0.1864 |

### Table S9. Sequence length distribution for unclassified sequences

| Statistic | Unclassified length (bp) |
|-----------|-------------------------|
| Minimum | 256 |
| Q1 (25th percentile) | 4,614 |
| Median | 12,516 |
| Q3 (75th percentile) | 112,671 |
| Maximum | 4,907,185 |

### Table S10. Forced Inc type distribution of unclassified sequences

| Inc type | Count | Percentage | Interpretation |
|----------|-------|------------|----------------|
| IncFII | 519 | 24.6% | Highly diverse group; boundary sequences expected |
| IncX1 | 318 | 15.1% | Compositionally heterogeneous; mosaic plasmids |
| IncN | 307 | 14.6% | Diverse backbone variants |
| IncI | 133 | 6.3% | Structural diversity in transfer regions |
| ColRNAI | 104 | 4.9% | Small plasmids; limited compositional signal |
| IncI1 | 86 | 4.1% | — |
| IncFIB | 79 | 3.7% | — |
| IncHI2 | 65 | 3.1% | — |
| IncI2 | 63 | 3.0% | — |
| IncX3 | 59 | 2.8% | — |
| IncHI1 | 56 | 2.7% | — |
| IncAC2 | 50 | 2.4% | — |
| IncA | 45 | 2.1% | — |
| Other (15 groups) | 225 | 10.7% | Combined remaining Inc/Rep groups |

### Table S11. Confidence score distribution for unclassified sequences

| Statistic | Confidence score |
|-----------|-----------------|
| Minimum | 20.0% |
| Q1 (25th percentile) | 36.6% |
| Median | 38.8% |
| Q3 (75th percentile) | 39.6% |
| Maximum | 39.9% |

Note: The majority of unclassified sequences (75%) have confidence scores between 36.6% and 39.9%, indicating they are borderline cases near the 40% threshold rather than fundamentally unclassifiable. Only a small minority have confidence scores below 30%.

---

## Biological interpretation

The 2,109 unclassified sequences share three distinguishing characteristics:

1. **Short length.** Median length of 12,516 bp versus 59,311 bp for classified sequences (4.7-fold shorter), consistent with small cryptic plasmids or assembly fragments that contain limited compositional signal for KNN classification.

2. **High divergence.** Median nearest-neighbour distance of 0.0206 versus 0.0067 for classified sequences (3.1-fold more distant), indicating these sequences occupy sparsely represented regions of tetranucleotide composition space.

3. **High individuality.** 1,991 unique pLIN codes among 2,109 sequences (94.4% unique), compared with 72.4% uniqueness for classified sequences, confirming that unclassified sequences are individually divergent rather than forming coherent unrecognised clusters.

These properties indicate that the unclassified fraction consists predominantly of short, individually divergent plasmids or fragments that lack the compositional complexity needed for reliable group assignment. Importantly, the novel Inc/Rep group discovery module (L6 analysis) identified 18 putative novel groups comprising 247 plasmids (11.7% of the unclassified fraction) from these low-confidence sequences, actively characterising the boundary of the classification system. The unclassifiable fraction is expected to decrease incrementally as new Inc/Rep groups are incorporated into the training set.

---

## Handling of unclassified sequences across methods

| Method | Behaviour for unclassifiable input | User notification |
|--------|-----------------------------------|-------------------|
| PlasmidFinder | No replicon hit returned; sequence silently excluded | None |
| pMLST | Ineligible if not in typed Inc group; no result returned | None |
| COPLA | "Unassigned" label; no further characterisation | Label only |
| MOB-typer | "Unknown" relaxase; reduced functional annotation | Label only |
| mge-cluster | Assigned to "noise" cluster; excluded from analysis | Cluster label |
| MOB-cluster | Always assigned a cluster; cluster ID may change between runs | None (no stability warning) |
| **pLIN** | **Forced assignment with explicit low-confidence flag, NN distance, and nearest reference ID** | **Confidence score, distance metric, nearest neighbour** |

pLIN provides the most informative handling of unclassifiable sequences among all compared methods: rather than silently excluding sequences or assigning unstable identifiers, it reports the best-match group alongside quantitative reliability metrics, enabling users to make informed decisions about downstream inclusion.
