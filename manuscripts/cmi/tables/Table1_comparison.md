# Table 1. Comparison of plasmid typing approaches for clinical infection control

| Feature | PlasmidFinder [5] | pMLST [6] | MOB-suite [7] | COPLA [8] | pLIN (this study) |
|---------|:-----------------:|:----------:|:--------------:|:---------:|:-----------------:|
| **Typing basis** | Replicon markers | Allelic profiles | Relaxase + replicon | Network clustering | 4-mer composition (LIN framework) |
| **Resolution levels** | 1 (Inc group) | 1 (sequence type) | 1 (MOB cluster) | 1 (PTU) | 6 (L1--L6) |
| **Discriminatory power (Simpson's D)** | 0.641^a^ | N/A | N/A | N/A | 0.985 |
| **Stable nomenclature across database updates** | Yes | Yes | No^b^ | No^b^ | Yes |
| **Integrated AMR profiling** | No | No | No | No | Yes (AMRFinderPlus) |
| **Lineage-level AMR risk profiles** | No | No | No | No | Yes |
| **Automated outbreak detection** | No | No | No | No | Yes (2-tier) |
| **Clinical risk stratification** | No | No | No | No | Yes (CRITICAL/HIGH/MODERATE) |
| **SNP sub-typing for outbreak confirmation** | No | No | No | No | Yes |
| **Multi-resolution investigation** | No | No | No | Partial | Yes (family to strain) |
| **Reference-free operation** | No | No | No | No | Yes |
| **Standard hardware (<30 min)** | Yes | Yes | Yes | Yes | Yes |

^a^ Calculated on the 8,077-plasmid training dataset (28 replicon groups) using Inc/rep group as the unit of classification.

^b^ MOB-suite and COPLA reassign cluster identifiers with each database release, preventing longitudinal comparison of codes across database versions.

N/A, not applicable or not reported.
