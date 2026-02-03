# pLIN: Plasmid Life Identification Number System

A hierarchical, permanent, reference-free classification system for bacterial plasmid genomes integrated with antimicrobial resistance gene surveillance.

## Overview

**pLIN** (plasmid Life Identification Number) is the first application of the Life Identification Number (LIN) framework to plasmid genomes. It assigns each plasmid a six-position hierarchical code (A.B.C.D.E.F) based on tetranucleotide composition distances and single-linkage clustering at six biologically calibrated thresholds.

### Key Features

- **Hierarchical**: 6 nested levels from family (~85% ANI) to strain (~99.9% ANI)
- **Permanent**: Codes never change when new plasmids are added
- **Reference-free**: No external database required — works on raw nucleotide sequences
- **AMR-integrated**: Full integration with NCBI AMRFinderPlus for resistance gene surveillance
- **ML-validated**: Nested cross-validation confirms compositional features predict Inc-group membership (F1=0.903)

## Dataset

- **6,346** complete plasmid genomes from NCBI RefSeq
- **3 Inc groups**: IncFII (n=4,581), IncN (n=1,064), IncX1 (n=701)
- **2,232** unique strain-level pLIN codes (Simpson's D = 0.979)

## Results Highlights

| Metric | Value |
|--------|-------|
| Total plasmids | 6,346 |
| Unique pLIN codes | 2,232 |
| Simpson's Diversity (D) | 0.979 |
| Inc-group concordance | 99.5% |
| ML best F1 (XGBoost) | 0.903 |
| AMR gene detections | 27,465 |
| Virulence detections | 5,834 |
| Plasmids with any hit | 84.2% |
| Carbapenemase detections | 1,490 |
| Colistin resistance (mcr) | 160 |

## Repository Structure

```
pLIN/
├── pLIN.ipynb                  # Main analytical notebook
├── pLIN_executed.ipynb         # Executed notebook with outputs
├── assign_pLIN.py              # pLIN assignment pipeline
├── run_amrfinder_all.sh        # AMRFinderPlus batch runner
├── integrate_pLIN_AMR.py       # pLIN + AMR integration script
├── generate_figures.py         # Publication figure generation
├── manuscript_complete.md      # Complete manuscript (Introduction, Methods, Results, Discussion)
├── manuscript_results.md       # Results section (standalone)
├── Data/
│   └── IncX_PLIN_thresholds_v0_python.yaml  # Threshold calibration data
├── output/
│   ├── pLIN_assignments.tsv                  # pLIN codes for all 6,346 plasmids
│   ├── amrfinder/
│   │   └── amrfinder_all_plasmids.tsv       # AMRFinderPlus raw detections
│   ├── integrated/
│   │   ├── pLIN_AMR_integrated.tsv          # Combined pLIN + AMR table
│   │   └── pLIN_lineage_AMR_summary.tsv     # Per-lineage AMR profiles
│   └── figures/
│       ├── Figure1_dataset_overview.png/pdf
│       ├── Figure2_AMR_prevalence.png/pdf
│       ├── Figure3_critical_AMR.png/pdf
│       ├── Figure4_pLIN_AMR_heatmap.png/pdf
│       ├── Figure5_AMR_burden_virulence.png/pdf
│       ├── Figure6_pLIN_hierarchy.png/pdf
│       └── Figure7_composite_manuscript.png/pdf
└── README.md
```

## Quick Start

### Requirements

```
Python >= 3.10
numpy
pandas
scipy
biopython
scikit-learn
xgboost
optuna
matplotlib
seaborn
```

### Install Dependencies

```bash
python -m venv .venv
source .venv/bin/activate
pip install numpy pandas scipy biopython scikit-learn xgboost optuna matplotlib seaborn
```

### Assign pLIN Codes

Place your FASTA files in `plasmid_sequences_for_training/<IncType>/fastas/` and run:

```bash
python assign_pLIN.py
```

Output: `output/pLIN_assignments.tsv`

### Run AMRFinderPlus Integration

Requires [AMRFinderPlus](https://github.com/ncbi/amr) installed:

```bash
# Run AMRFinderPlus on all plasmids
bash run_amrfinder_all.sh

# Integrate with pLIN codes
python integrate_pLIN_AMR.py
```

### Generate Figures

```bash
python generate_figures.py
```

## pLIN Code Format

Each plasmid receives a six-position code: `A.B.C.D.E.F`

| Position | Level | Distance Threshold | ANI Equivalent |
|----------|-------|-------------------|----------------|
| A | Family | d ≤ 0.150 | ~85% |
| B | Subfamily | d ≤ 0.100 | ~90% |
| C | Cluster | d ≤ 0.050 | ~95% |
| D | Subcluster | d ≤ 0.020 | ~98% |
| E | Clone complex | d ≤ 0.010 | ~99% |
| F | Strain | d ≤ 0.001 | ~99.9% |

Example: pLIN `1.1.1.7.30.567` = Family 1, Subfamily 1, Cluster 1, Subcluster 7, Clone 30, Strain 567

## Plasmid Sequences

FASTA sequences (695 MB total) are not included in this repository due to size. All plasmid sequences were obtained from NCBI RefSeq and can be downloaded using the accession numbers listed in `output/pLIN_assignments.tsv`.

## Citation

If you use pLIN in your research, please cite:

> Xavier, B. (2025). pLIN: A Plasmid Life Identification Number System for Hierarchical, Permanent Classification of Bacterial Plasmids Integrated with Antimicrobial Resistance Gene Surveillance.

## License

MIT License
