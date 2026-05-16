# Table 2: pLIN hierarchical threshold definitions

Six cosine distance thresholds define the hierarchical levels of the pLIN classification system. Thresholds were initially calibrated using the nearest-neighbour distance distribution of 178 IncX-like reference plasmids and validated against FastANI across 4,970 plasmid pairs from all 20 Gram-negative Inc groups.

| Level | Cosine distance (d) | Approximate ANI | Taxonomic analogy | Clinical interpretation | Calibration quantile |
|-------|:---:|:---:|---|---|---|
| L1 | <=0.150 | ~85% | Plasmid superfamily | Broad plasmid family identification; incompatibility group-level | 99th percentile |
| L2 | <=0.100 | ~90% | Major lineage | Epidemiological lineage grouping; regional surveillance | 95th--99th percentile |
| L3 | <=0.050 | ~95% | Species-level cluster | Regional surveillance comparisons; epidemic context | 95th percentile (d=0.048) |
| L4 | <=0.020 | ~98% | Sublineage | Inter-hospital transmission tracking | 75th percentile (d=0.025) |
| L5 | <=0.010 | ~99% | Clone complex | Intra-hospital transmission investigation; ward-level spread | Median (d=0.011) |
| L6 | <=0.001 | ~99.9% | Strain / outbreak level | Outbreak confirmation; infection control action trigger | 25th percentile (d=0.001) |

**FastANI validation (n=4,970 pairs):** At L6 (d<=0.001), median ANI = 99.9% (IQR 98.8--100.0%). Overall Spearman rho = -0.348 (P < 10^-141). Significant per-group correlations (P < 0.001) were observed for 15 of 20 Gram-negative Inc groups; the five non-significant groups (IncX1, IncFII, IncFIBK, IncX4, IncX3) exhibited narrow within-group ANI ranges (>89%), producing a floor effect consistent with high compositional homogeneity.
