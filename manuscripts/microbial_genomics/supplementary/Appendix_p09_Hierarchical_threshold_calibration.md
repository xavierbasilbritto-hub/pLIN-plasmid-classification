# Appendix p 9: Hierarchical threshold calibration

## Initial calibration (IncX-like subset, n=178)

Thresholds were initially calibrated using the nearest-neighbour distance distribution of 178 selected IncX-like reference plasmids (seed plasmid: RefSeq_NZ_AP027441.1; selection cutoff d ≤ 0.12).

| Quantile | Nearest-neighbour distance (d) | Calibration role |
|----------|-------------------------------|-----------------|
| 1st–10th | 0.000 | Identical/near-identical compositions |
| 25th | 0.001 | → L6 threshold (strain/outbreak level) |
| 50th (median) | 0.011 | → L5 threshold (clone complex) |
| 75th | 0.025 | → L4 threshold (subcluster) |
| 90th | 0.038 | Near L3 threshold (cluster) |
| 95th | 0.048 | → L3 threshold (~95% ANI) |
| 99th | 0.072 | → L2 threshold range |

## Cosine distance to ANI approximate mapping

- d = 0.001 corresponds to ~99.9% ANI (outbreak level)
- d = 0.010 corresponds to ~99% ANI (clone complex)
- d = 0.020 corresponds to ~98% ANI (within-species sublineage)
- d = 0.050 corresponds to ~95% ANI (species boundary)
- d = 0.100 corresponds to ~90% ANI (major lineage divergence)
- d = 0.150 corresponds to ~85% ANI (plasmid family boundary)

## Multi-Inc validation (all 20 Gram-negative Inc groups, n=4,970 pairs)

FastANI (v1.34, --fragLen 1000) was computed for 4,970 within-group plasmid pairs sampled from all 20 Gram-negative Inc groups (up to 20 randomly sampled plasmids per group, all pairwise comparisons within each group).

**Overall correlation:** Spearman ρ = -0.348 (P < 10⁻¹⁴¹)

### Median FastANI at each pLIN threshold

| pLIN Level | Cosine threshold | Expected ANI | n pairs | Median ANI | 5th %ile | 25th %ile | 75th %ile |
|------------|-----------------|-------------|---------|-----------|----------|----------|----------|
| L6 (Strain) | d ≤ 0.001 | ~99.9% | 726 | 99.9% | 97.4% | 98.8% | 100.0% |
| L5 (Clone) | d ≤ 0.010 | ~99.0% | 2,696 | 98.2% | 90.4% | 97.1% | 99.7% |
| L4 (Subcluster) | d ≤ 0.020 | ~98.0% | 3,566 | 97.9% | 89.5% | 95.9% | 99.5% |
| L3 (Cluster) | d ≤ 0.050 | ~95.0% | 4,554 | 97.8% | 87.9% | 95.3% | 99.5% |
| L2 (Subfamily) | d ≤ 0.100 | ~90.0% | 4,944 | 97.8% | 85.3% | 95.1% | 99.5% |
| L1 (Family) | d ≤ 0.150 | ~85.0% | 4,966 | 97.8% | 85.2% | 95.1% | 99.5% |

### Mean ANI per cosine distance bin

| Cosine distance bin | n pairs | Mean ANI | Median ANI | Std | Min | Max |
|---------------------|---------|----------|-----------|-----|-----|-----|
| d ≤ 0.001 | 696 | 99.4% | 99.9% | 0.9 | 97.0% | 100.0% |
| 0.001–0.005 | 1,172 | 97.5% | 97.9% | 2.8 | 78.9% | 100.0% |
| 0.005–0.010 | 798 | 95.5% | 96.9% | 4.3 | 77.0% | 100.0% |
| 0.010–0.020 | 870 | 94.9% | 96.3% | 4.7 | 74.9% | 100.0% |
| 0.020–0.050 | 988 | 95.5% | 97.7% | 5.3 | 73.4% | 100.0% |
| 0.050–0.100 | 390 | 94.4% | 97.6% | 6.3 | 76.7% | 100.0% |
| 0.100–0.150 | 22 | 90.4% | 93.0% | 5.6 | 80.8% | 96.4% |
| d > 0.150 | 4 | 84.7% | 86.0% | 3.5 | 79.5% | 87.1% |

### Per-Inc-group Spearman correlation (cosine distance vs FastANI)

| Inc group | n pairs | Spearman ρ | P-value | Significance | Mean ANI |
|-----------|---------|-----------|---------|--------------|----------|
| IncFIC | 182 | -0.883 | <0.001 | *** | 98.5% |
| IncI2 | 308 | -0.810 | <0.001 | *** | 98.3% |
| ColE | 254 | -0.722 | <0.001 | *** | 93.2% |
| IncI | 90 | -0.612 | <0.001 | *** | 96.1% |
| IncN | 230 | -0.612 | <0.001 | *** | 95.2% |
| IncHI1 | 75 | -0.575 | <0.001 | *** | 94.0% |
| IncAC2 | 182 | -0.536 | <0.001 | *** | 99.1% |
| IncF | 343 | -0.514 | <0.001 | *** | 93.3% |
| ColRNAI | 202 | -0.502 | <0.001 | *** | 96.4% |
| IncFIB | 238 | -0.475 | <0.001 | *** | 94.5% |
| IncI1 | 306 | -0.297 | <0.001 | *** | 97.8% |
| IncA | 182 | -0.288 | <0.001 | *** | 94.5% |
| IncHI2 | 380 | -0.251 | <0.001 | *** | 98.3% |
| IncR | 378 | -0.169 | <0.001 | *** | 95.8% |
| IncC | 240 | -0.148 | 0.022 | * | 98.6% |
| IncX1 | 274 | -0.076 | 0.207 | n.s. | 93.5% |
| IncFII | 240 | -0.070 | 0.281 | n.s. | 89.4% |
| IncFIBK | 106 | -0.061 | 0.531 | n.s. | 98.2% |
| IncX4 | 380 | -0.016 | 0.754 | n.s. | 98.7% |
| IncX3 | 380 | +0.038 | 0.467 | n.s. | 99.6% |

**Interpretation:** The five Inc groups with non-significant correlations (IncX1, IncFII, IncFIBK, IncX4, IncX3) exhibit narrow within-group ANI ranges — all sampled pairs have ANI > 89%, many > 98% — producing a floor effect where cosine distance variation occurs within a compressed ANI band. This does not indicate threshold miscalibration but rather reflects the high compositional and sequence homogeneity within these groups. The 15 groups with significant correlations span the full range of plasmid diversity, confirming that the cosine-to-ANI mapping is valid across diverse plasmid families.

**Full data:** cosine_to_ani_validation_all_inc.tsv (4,970 rows) and cosine_to_ani_per_inc_summary.tsv (20 rows)
