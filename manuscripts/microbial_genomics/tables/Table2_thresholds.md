# Table 2. pLIN hierarchical threshold definitions and clustering results for the 8,077-sample (8,056 unique) training dataset.

| Level | Designation | Cosine distance threshold (d) | Approximate ANI equivalent | Clusters (n) | Median cluster size | Maximum cluster size | Singletons (n) |
|-------|-------------|-------------------------------|---------------------------|---------------|---------------------|----------------------|-----------------|
| L1 | Family | <= 0.150 | ~85% | 1 | 6,998 | 6,998 | 0 |
| L2 | Subfamily | <= 0.100 | ~90% | 1 | 6,998 | 6,998 | 0 |
| L3 | Cluster | <= 0.050 | ~95% | 7 | 1 | 6,990 | 5 |
| L4 | Subcluster | <= 0.020 | ~98% | 38 | 1 | 6,947 | 27 |
| L5 | Clone complex | <= 0.010 | ~99% | 117 | 1 | 6,737 | 82 |
| L6 | Strain | <= 0.001 | ~99.9% | 3,073 | 1 | 869 | 2,335 |

Thresholds were calibrated against the nearest-neighbour distance distribution of 178 IncX-like reference plasmids (seed: RefSeq NZ_AP027441.1). Key calibration quantiles: 50th percentile d = 0.011, 90th percentile d = 0.038, 99th percentile d = 0.072. ANI equivalents were validated using FastANI v1.34 on 4,970 within-group plasmid pairs from the 20 Gram-negative Inc groups (Spearman rho = -0.348, P < 10^-141; 15/20 groups significant at P < 0.05). No FastANI validation was performed for the Gram-positive, *Acinetobacter*, or *Pseudomonas* groups. At d <= 0.001, the median FastANI was 99.9%.
