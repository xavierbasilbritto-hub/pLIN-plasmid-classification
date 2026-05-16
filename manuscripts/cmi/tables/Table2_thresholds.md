# Table 2. pLIN hierarchical thresholds, approximate ANI equivalence, and clinical interpretation

| Level | Cosine distance threshold (d) | Approximate ANI | Median ANI (n=4,970 pairs) | Clinical interpretation | Infection control use case |
|:-----:|:----------------------------:|:---------------:|:--------------------------:|------------------------|---------------------------|
| L1 | <=0.150 | ~85% | 97.8% | Plasmid superfamily | Broad plasmid family identification; epidemiological context |
| L2 | <=0.100 | ~90% | 97.8% | Major lineage | National/international lineage-level surveillance |
| L3 | <=0.050 | ~95% | 97.8% | Species-level cluster | Regional surveillance; comparison across healthcare networks |
| L4 | <=0.020 | ~98% | 97.9% | Sublineage | Inter-hospital transmission tracking; referral pathway analysis |
| L5 | <=0.010 | ~99% | 98.2% | Clone complex | Intra-hospital transmission investigation; ward-level spread |
| L6 | <=0.001 | ~99.9% | 99.9% | Strain / outbreak level | Outbreak confirmation; direct transmission evidence |

Thresholds were initially calibrated on the nearest-neighbour distance distribution of 178 IncX-like plasmids and validated across the 20 Gram-negative Inc groups using FastANI (v1.34, --fragLen 1000) on 4,970 within-group plasmid pairs (Spearman rho = -0.348, P < 10^-141; 15/20 groups significant). FastANI validation was not performed for Gram-positive, *Acinetobacter*, or *Pseudomonas* groups. The pLIN code format at each level is cumulative (e.g., 1.1.2.15.48.671 represents L1.L2.L3.L4.L5.L6).
