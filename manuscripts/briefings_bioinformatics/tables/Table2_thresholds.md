# Table 2. pLIN hierarchical threshold definitions, ANI calibration, and clustering results at each level.

## A. Threshold definitions and biological interpretation

| Level | Bin | Cosine distance threshold (d) | ANI equivalent | Biological interpretation | Typical application |
|-------|-----|-------------------------------|----------------|---------------------------|---------------------|
| L1 | A | <= 0.150 | ~85% | Family | Broad evolutionary comparisons |
| L2 | B | <= 0.100 | ~90% | Subfamily | Plasmid family substructure |
| L3 | C | <= 0.050 | ~95% | Cluster | Lineage tracking |
| L4 | D | <= 0.020 | ~98% | Subcluster | Lineage-level surveillance |
| L5 | E | <= 0.010 | ~99% | Clone complex | Clonal expansion detection |
| L6 | F | <= 0.001 | ~99.9% | Strain / Outbreak | Outbreak investigation |

## B. Clustering results at each threshold for the full training dataset (n = 8,077 samples, 28 groups)

| Bin | Level | Distance threshold (d) | Clusters (n) | Median cluster size | Max cluster size | Singletons (n) | Singletons (%) |
|-----|-------|------------------------|---------------|---------------------|------------------|-----------------|----------------|
| A | Family (L1) | <= 0.150 | 1 | 8,077 | 8,077 | 0 | 0.0% |
| B | Subfamily (L2) | <= 0.100 | 2 | -- | -- | 0 | 0.0% |
| C | Cluster (L3) | <= 0.050 | 9 | 1 | -- | -- | -- |
| D | Subcluster (L4) | <= 0.020 | 44 | 1 | -- | -- | -- |
| E | Clone complex (L5) | <= 0.010 | 131 | 1 | -- | -- | -- |
| F | Strain (L6) | <= 0.001 | 3,073 | 1 | 869 | 2,335 | 76.0% |

## C. Comparison: training-only vs. full reference database (n = 79,305)

| Bin | Level | Training only (n = 8,077) | Full database (n = 77,196^a^) | Fold increase |
|-----|-------|---------------------------|-------------------------------|---------------|
| A | Family (L1) | 1 | 33 | 33x |
| B | Subfamily (L2) | 2 | 86 | 43x |
| C | Cluster (L3) | 9 | 447 | 50x |
| D | Subcluster (L4) | 44 | 2,772 | 63x |
| E | Clone complex (L5) | 131 | 5,540 | 42x |
| F | Strain (L6) | 3,073 | 57,886 | 18.8x |

^a^ 77,196 classified sequences (79,305 total minus 2,109 flagged as Unknown/Novel at <40% confidence).

## D. FastANI validation of cosine-to-ANI mapping (n = 4,970 pairs)

| Metric | Value |
|--------|-------|
| Total pairwise comparisons | 4,970 |
| Overall Spearman rho | -0.348 |
| P-value | < 10^-141 |
| Inc groups with significant correlation (P < 0.05) | 15 / 20 |
| Inc groups with strong correlation (rho < -0.2) | 13 / 20 |
| Median ANI at d <= 0.001 | 99.9% |
| Strongest per-group rho | IncFIC (-0.88), IncI2 (-0.81), ColE (-0.72), IncN (-0.61) |

**Notes:**

- Thresholds were initially calibrated against the nearest-neighbour distance distribution of 178 IncX-like reference plasmids (seed plasmid: RefSeq_NZ_AP027441.1, distance cutoff d <= 0.12). Nearest-neighbour quantiles: 50th percentile d = 0.011, 90th percentile d = 0.038, 99th percentile d = 0.072.
- FastANI validation used --fragLen 1000 on up to 20 plasmids per Inc group (all pairwise comparisons within each group).
- Five groups (IncX3, IncX4, IncFIBK, IncFII, IncX1) showed non-significant cosine-ANI correlations due to range restriction (all within-group ANI values > 93%), producing a floor effect.

**Abbreviations:** ANI, average nucleotide identity; d, cosine distance.
