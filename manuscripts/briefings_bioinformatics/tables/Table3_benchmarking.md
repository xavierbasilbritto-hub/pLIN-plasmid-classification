# Table 3. Comprehensive benchmarking results for the pLIN framework.

## A. Discriminatory power (Simpson's Index of Diversity)

| Classification system | D value | Unique codes (n) | Dataset |
|-----------------------|---------|-------------------|---------|
| pLIN at L6 (strain) | **0.985** | 3,073 | 8,077 samples, 28 groups |
| pLIN at L5 (clone complex) | 0.963 | 131 | 8,077 samples, 28 groups |
| pLIN at L4 (subcluster) | 0.878 | 44 | 8,077 samples, 28 groups |
| Inc/rep typing alone | 0.641 | 28 | 8,077 samples, 28 groups |
| Improvement (pLIN L6 / Inc) | **1.54-fold** | -- | -- |

### Per-Inc-group discriminatory power at L6

| Inc group | n | D value | Unique pLIN codes |
|-----------|---|---------|-------------------|
| IncX1 | 705 | 0.995 | 420 |
| IncN | 1,097 | 0.976 | 431 |
| IncFII | 4,629 | 0.962 | 1,421 |

## B. Group concordance at strain level (L6)

| Metric | Value |
|--------|-------|
| Total unique pLIN codes | 3,073 |
| Single-group codes | 3,006 (97.8%) |
| Multi-group codes | 67 (2.2%) |
| Largest mixed-group cluster | pLIN 1327 (n = 869, 6 Inc groups) |
| Singletons | 2,335 (76.0%) |

## C. Machine learning cross-validation (nested CV: 3 outer folds, 2 inner folds, Optuna TPE, 10 trials/fold, Hyperband pruning)

| Model | Mean weighted F1 | SD | Min | Max |
|-------|-------------------|----|-----|-----|
| XGBoost | **0.896** | 0.009 | 0.887 | 0.905 |
| Gradient Boosting | 0.893 | 0.007 | 0.886 | 0.900 |
| Random Forest | 0.874 | 0.019 | 0.855 | 0.893 |
| Logistic Regression | 0.866 | 0.010 | 0.856 | 0.876 |

### Top 5 discriminative features (consensus across tree-based models)

| Rank | Feature | Mean importance | Biological relevance |
|------|---------|-----------------|----------------------|
| 1 | tri_TAG | 0.114 | Amber stop codon; codon usage signal |
| 2 | tri_GCG | 0.099 | CpG island proxy; methylation signature |
| 3 | at_content | 0.076 | Overall base composition |
| 4 | gc_content | 0.053 | GC%; backbone vs. accessory content |
| 5 | di_AA | 0.050 | Poly-A tracts; regulatory elements |

## D. KNN Inc-group classifier performance

| Metric | Value |
|--------|-------|
| Algorithm | K-nearest neighbours |
| k | 5 |
| Distance metric | Cosine |
| Weighting | Distance-weighted |
| Feature space | 256 tetranucleotide frequencies |
| Training samples | 8,077 (8,056 unique) |
| Groups | 28 |
| Cross-validation | 5-fold |
| Overall accuracy | **91.1%** |
| Classification rate (reference DB) | 97.3% (69,140 / 71,249) |
| Unknown/Novel (< 40% confidence) | 2,109 (2.7%) |

### Confidence-stratified accuracy

| Confidence threshold | % of reference sequences | Estimated accuracy |
|----------------------|--------------------------|--------------------|
| >= 95% | 43.3% | 96.6% |
| >= 80% | 57.0% | >= 86.3% |
| < 40% (excluded) | 3.3% | -- |

## E. Scalability benchmarks (Apple M-series laptop, 16--32 GB RAM)

| Metric | Value |
|--------|-------|
| Total sequences processed | 79,305 |
| Total runtime | < 30 minutes |
| Phase 1 (4-mer vectors, 71,249 seqs) | 23.5 min (~52 seq/sec) |
| Phase 2 (KNN classification) | 7 sec (~10,272 seq/sec) |
| Phase 3 (per-group clustering, 28 groups) | 4.2 min |
| Phase 4 (output) | < 1 sec |
| Maximum group memory (IncFII, 34,036 seqs) | **2.3 GB** |
| Global clustering memory (avoided) | ~24 GB |
| Memory reduction factor | **10.4x** |
| Unique pLIN codes resolved | 57,886 |
| Fold increase over training only | 13.5x |

## F. Outbreak cross-validation (74 plasmids, 27 published studies, 13 countries)

| Validation criterion | Result |
|---------------------|--------|
| Total outbreak plasmids tested | 74 |
| Independent published studies | 27 |
| Countries represented | 13 (4 continents) |
| Resistance mechanisms tested | 7 (KPC, NDM, OXA-48, VIM, IMP, mcr, CTX-M) |
| Unique pLIN codes assigned | 42 |
| High-confidence classifications (>=60%) | 63/74 (85.1%) |
| Known high-risk lineage matches | 3 (pLIN 671 [KPC-2], pLIN 860 [MDR hub], pLIN 1688 [OXA-48]) |
| Intra-study outbreak clusters | 9 |
| Cross-continent code matches | Yes (OXA-48: 3 countries; mcr-1: 3 countries) |
| Mean nearest-neighbour distance | 0.0062 |

### Individual study results

| Study | Accession | Assigned pLIN | Inc type | Confidence | NN distance | Match |
|-------|-----------|---------------|----------|------------|-------------|-------|
| Yao 2023 [17] | CP104944 | 1.1.2.15.48.671 | IncN | 100% | 0.0000 | Known KPC-2 lineage |
| Weber 2019 [18] | 12 plasmids | 9 unique L6, 1 L3 cluster | IncN/IncAC2/IncC/IncFII | 41--100% | 0.0000--0.0039 | Correct hierarchical grouping |
| Marimuthu 2022 [19] | MN542377 | 1.1.2.15.48.2455 | IncN | 60.9% | 0.0058 | Novel variant, correct neighbourhood |
| Roberts 2020 [20] | CP022533 | 1.1.2.15.48.860 | IncHI2 | 100% | 0.0004 | Known MDR hub lineage |

## G. AMR integration summary (AMRFinderPlus v4.2.5, n = 6,998 Gram-negative plasmids, 20 Inc groups)

| Metric | Value |
|--------|-------|
| Total gene detections | 64,891 |
| AMR gene detections | 29,583 |
| Virulence factor detections | 6,286 |
| Stress response detections | 29,022 |
| Plasmids with any hit | 5,816 / 6,998 (83.1%) |
| Plasmids with AMR genes | 4,657 / 6,998 (66.5%) |
| Carbapenemase detections | 1,635 |
| ESBL detections | 1,804 |
| Colistin resistance (mcr) detections | 204 |
| PMQR detections | 2,315 |

### High-risk pLIN lineages

| pLIN code | Inc type(s) | n | Mean AMR genes | Key determinants |
|-----------|-------------|---|----------------|------------------|
| 671 | IncN | 90 | 13.2 | blaKPC-2 (100%), blaTEM-1 (100%), aph(3'')-Ib (98%) |
| 860 | IncN/IncHI2/IncFII/IncHI1/IncX1 | 142 | 14.4 | sul1 (61%), floR (56%), mcr-1.1 (37%), 44.4% mcr total |
| 1432 | IncFII | 21 | 17.1 | blaNDM-5 (71%), blaTEM-1 (95%), sul1 (100%) |

**Abbreviations:** AMR, antimicrobial resistance; CV, cross-validation; D, Simpson's Index of Diversity; DB, database; ESBL, extended-spectrum beta-lactamase; KNN, K-nearest neighbours; MDR, multidrug-resistant; NN, nearest neighbour; PMQR, plasmid-mediated quinolone resistance; SD, standard deviation; TPE, tree-structured Parzen estimator.
