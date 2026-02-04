# pLIN: Plasmid Life Identification Number System

A hierarchical, permanent, reference-free classification system for bacterial plasmid genomes integrated with antimicrobial resistance gene surveillance.

## Overview

**pLIN** (plasmid Life Identification Number) is the first application of the Life Identification Number (LIN) framework to plasmid genomes. It assigns each plasmid a six-position hierarchical code (A.B.C.D.E.F) based on tetranucleotide composition distances and single-linkage clustering at six biologically calibrated thresholds.

### Key Features

- **Hierarchical**: 6 nested levels from family (~85% ANI) to strain (~99.9% ANI)
- **Permanent**: Codes never change when new plasmids are added
- **Reference-free**: No external database required — works on raw nucleotide sequences
- **AMR-integrated**: Full integration with NCBI AMRFinderPlus for resistance gene surveillance
- **Inc Auto-Detection**: KNN classifier (96.1% accuracy) identifies Inc groups from k-mer composition
- **Mobility Prediction**: Classifies plasmids as conjugative, mobilizable, or non-mobilizable
- **Outbreak Detection**: Flags clonal clusters sharing identical pLIN codes and AMR profiles
- **Adaptive Thresholds**: Per-Inc-group distance threshold calibration from training data
- **Multi-Linkage**: Selectable clustering methods (single, complete, average, weighted)
- **Interactive GUI**: Streamlit web app with 6 analysis tabs — no command line needed
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
├── plin_app.py                 # Streamlit GUI (main interactive app)
├── assign_pLIN.py              # Batch pLIN assignment pipeline
├── build_inc_centroids.py      # Train Inc group KNN classifier
├── integrate_pLIN_AMR.py       # pLIN + AMR integration script
├── generate_figures.py         # Publication figure generation
├── create_architecture_pptx.py # PowerPoint architecture generator
├── test_pLIN.py                # Test pLIN on 22 plasmids
├── test_cladogram.py           # Test cladogram generation
├── test_integrate_and_cladogram.py  # Test AMR + cladogram
├── train_nt_classifier.py      # Train Nucleotide Transformer probes
├── launch_pLIN.bat             # Windows one-click launcher
├── launch_pLIN.command         # macOS one-click launcher
├── launch_pLIN.sh              # Linux one-click launcher
├── setup.sh / setup.bat        # One-command setup
├── run_all.sh / run_all.bat    # Run complete pipeline
├── run_amrfinder_all.sh        # AMRFinderPlus batch runner
├── requirements.txt            # Python dependencies
├── pLIN.ipynb                  # Main analytical notebook
├── manuscript_complete.md      # Complete manuscript
├── LICENSE                     # GPL-3.0 + Citation clause
├── CITATION.cff               # GitHub citation metadata
├── Data/
│   ├── inc_classifier.npz      # KNN classifier model (4.3 MB)
│   ├── inc_centroids.npz       # Inc group centroid profiles
│   └── IncX_PLIN_thresholds_v0_python.yaml
├── test_plasmids/              # Test FASTA files (IncX, IncFII, IncH)
├── output/
│   ├── pLIN_Tool_Architecture.pptx  # 11-slide architecture presentation
│   ├── pLIN_assignments.tsv
│   ├── integrated/
│   ├── figures/
│   └── test/                   # Test output (cladograms, AMR results)
└── README.md
```

## Quick Start

### One-Click Launch (Recommended)

Download the repository and double-click the launcher for your platform:

| Platform | Launcher File | How to Run |
|----------|--------------|------------|
| **Windows** | `launch_pLIN.bat` | Double-click the file |
| **macOS** | `launch_pLIN.command` | Double-click the file |
| **Linux** | `launch_pLIN.sh` | Run `chmod +x launch_pLIN.sh && ./launch_pLIN.sh` |

The launcher will automatically:
1. Check for Python 3.10+
2. Create a virtual environment
3. Install all dependencies
4. Launch the pLIN web app in your browser

**No command-line experience required** — just download and double-click.

### Manual Setup

```bash
# 1. Clone and setup
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
cd pLIN-plasmid-classification
pip install -r requirements.txt

# 2. Launch the web app
streamlit run plin_app.py
```

Then open your browser to `http://localhost:8501`, upload FASTA files, and click **Run Analysis**.

The GUI provides:
- **Overview** — pLIN system description, threshold table, analysis parameters
- **Results** — Interactive table with pLIN codes, Inc group, mobility, AMR data
- **Cladogram** — 4 visualization types (rectangular, circular, heatmap, AMR-annotated)
- **AMR Analysis** — Gene prevalence, drug class breakdown, critical gene alerts
- **Epidemiology** — Mobility prediction, outbreak detection, dissemination risk
- **Export** — Download TSV, PNG, PDF, ZIP bundle

### Command-Line Pipeline

```bash
# Setup
bash setup.sh                      # Creates venv, installs deps

# Add FASTA files to plasmid_sequences_for_training/{IncFII,IncN,IncX1}/fastas/

# Run full pipeline
bash run_all.sh

# Or individual steps:
source .venv/bin/activate
python assign_pLIN.py              # Assign pLIN codes
bash run_amrfinder_all.sh          # Run AMRFinderPlus
python integrate_pLIN_AMR.py       # Integrate pLIN + AMR
python generate_figures.py         # Generate figures
```

### Requirements

- Python >= 3.10
- AMRFinderPlus (optional, for AMR gene detection):
  ```bash
  conda install -c bioconda -c conda-forge ncbi-amrfinderplus
  amrfinder -u  # update database
  ```

Python packages (installed automatically):
```
numpy, pandas, scipy, biopython, scikit-learn, matplotlib, seaborn, streamlit, plotly
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

## Advanced Features

### Inc Group Auto-Detection
A KNN classifier (k=5, cosine distance, distance-weighted) trained on 6,346 plasmids predicts Inc group membership from 4-mer composition with 96.1% cross-validation accuracy. No replicon gene BLAST required.

### Adaptive Thresholds
Calibrates pLIN distance thresholds per Inc group using quantile-based analysis of within-group distance distributions. Addresses the limitation that fixed thresholds may not optimally separate lineages across diverse Inc groups.

### Mobility Prediction
Scans AMRFinderPlus output for conjugation (tra/trb) and mobilization (mob) gene markers to classify each plasmid as:
- **Conjugative** — has transfer genes, can self-transfer (highest risk)
- **Mobilizable** — has mob genes, needs helper plasmid
- **Non-mobilizable** — no detectable transfer machinery

### Outbreak Detection
Automatically flags groups of plasmids sharing identical pLIN strain codes (F-level) AND the same AMR resistance profile, indicating potential clonal spread. Risk levels: HIGH (3+ shared AMR genes) / MODERATE (1-2).

### Multi-Linkage Clustering
Supports single (default), complete, average, and weighted linkage methods. Complete linkage produces tighter clusters and reduces the chaining artifact inherent to single linkage.

### Nucleotide Transformer LLM (Optional)
Optional integration with [InstaDeep's Nucleotide Transformer](https://github.com/instadeepai/nucleotide-transformer), a genomic foundation model pre-trained on DNA sequences. When enabled, the NT model provides an independent LLM-based prediction of Inc group and AMR drug classes alongside the traditional KNN classifier.

**How it works:**
1. Plasmid sequences are split into overlapping 5,000 bp chunks
2. Each chunk is embedded by the NT model (6-mer tokenization, transformer encoding)
3. Chunk embeddings are mean-pooled to produce a single vector per plasmid
4. Linear probes (trained on your data) predict Inc group and AMR classes

**Setup:**
```bash
# Install optional dependencies
pip install transformers torch

# Train probes from your labeled training data (one-time)
python train_nt_classifier.py --model 50m
```

Then enable "Use Nucleotide Transformer (LLM)" in the GUI. Model variants: 50M (fast), 100M, 250M, 500M parameters. CPU and GPU (CUDA/MPS) supported.

## Plasmid Sequences

FASTA sequences (695 MB total) are not included in this repository due to size. All plasmid sequences were obtained from NCBI RefSeq and can be downloaded using the accession numbers listed in `output/pLIN_assignments.tsv`.

## Citation (MANDATORY)

**Any use of this software, its algorithms, outputs, or derivative works in publications, presentations, reports, theses, or other academic/commercial work REQUIRES citation.** This is a binding condition of the license.

> Xavier, B. (2025). **pLIN: A Plasmid Life Identification Number System for Hierarchical, Permanent Classification of Bacterial Plasmids Integrated with Antimicrobial Resistance Gene Surveillance.** GitHub: https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification

### BibTeX

```bibtex
@software{xavier2025plin,
  author       = {Britto, Basil Xavier},
  title        = {{pLIN: A Plasmid Life Identification Number System for
                   Hierarchical, Permanent Classification of Bacterial
                   Plasmids Integrated with Antimicrobial Resistance Gene
                   Surveillance}},
  year         = {2025},
  url          = {https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification},
  version      = {1.0.0},
}
```

### Citation Requirements

1. **Publications**: Any journal article, conference paper, preprint, thesis, dissertation, poster, or presentation that uses pLIN MUST include the citation above.
2. **Derivative works**: Any software, tool, or pipeline incorporating pLIN code or algorithms MUST include attribution: *"Based on pLIN by Basil Xavier Britto (2025)"* with a link to this repository.
3. **Generated data**: Any publicly distributed dataset or results produced using pLIN MUST acknowledge the tool and cite the reference.
4. **No exception**: Failure to cite constitutes a violation of the license terms.

## Copyright

Copyright (C) 2025 Basil Xavier Britto. All rights reserved.

## License

This project is licensed under the **GNU General Public License v3.0** with an **additional mandatory citation clause** under Section 7 of the GPL.

This means:
- You **may** use, modify, and redistribute this software
- You **must** cite the original work in any publication or derivative work
- Any modified version **must** also be open-source under GPL-3.0
- The citation requirement and copyright notices **must not** be removed

See [LICENSE](LICENSE) for full terms.
