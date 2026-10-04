# pLIN: Plasmid Lineage Identification Number

Software: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23138103.svg)](https://doi.org/10.5281/zenodo.23138103) Database: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23126057.svg)](https://doi.org/10.5281/zenodo.23126057)

A permanent, multi-level nomenclature for bacterial plasmids, from backbone families to outbreak clones, with integrated antimicrobial resistance (AMR) gene surveillance.

**Author:** Basil Britto Xavier | **License:** GPL-3.0-or-later | **How to cite:** see [CITATION.cff](CITATION.cff)

---

## Overview

pLIN gives every plasmid a **six-level code** (for example `169.178.183.208.209.2438`) that never changes as the database grows. The levels run from broad to fine:

| Level | Meaning | How it is assigned |
|---|---|---|
| L1 | Backbone family | shared protein families (containment >= 0.40) |
| L2 | Backbone group | shared protein families (containment >= 0.60) |
| L3 | Shared backbone | shared k-mers (containment >= 0.50) |
| L4 | Backbone variant | shared k-mers (containment >= 0.80) |
| L5 | Lineage | symmetric k-mer similarity >= 0.80 to its nearest relative |
| L6 | Near-identical (outbreak clone) | symmetric k-mer similarity >= 0.95 to its nearest relative |

A plasmid with a near relative in the database (L5 similarity) copies that relative's code, so near-identical plasmids, mutants and outbreak isolates stay together. A plasmid without one is placed in the backbone levels by a founder rule and starts a new lineage. Existing codes are never renumbered, and a plasmid already in the database always gets back its published code. The scheme (pLIN v4.1) was evaluated in a pre-registered study; see [Validation](#validation).

### Key Features

| Category | Features |
|----------|----------|
| **Classification** | Permanent six-level pLIN codes (L1-L6); replicon (Inc/Rep) group prediction for 29 groups; flagging of new (provisional) codes |
| **AMR Surveillance** | AMRFinderPlus integration (AMR + stress + virulence genes), critical gene alerts, drug class analysis |
| **Genomic Analysis** | Mash/MinHash ANI estimation, FastANI, minimap2 SNP sub-typing within lineages |
| **Epidemiology** | Plasmid mobility prediction (MOB-suite + AMRFinderPlus), outbreak detection, temporal outbreak clustering (30-day window) |
| **Host Inference** | CRISPR spacer-based host prediction (MinCED + BLAST+) |
| **Visualization** | Interactive Streamlit GUI, cladograms, heatmaps, Plotly charts |
| **Deployment** | Cross-platform (macOS/Windows/Linux), Docker support, one-click launchers |

### Validation

Pre-registered on GitHub before the test data were drawn ([v4.1 pre-registration](../../blob/preregistration-v4/output/backbone_v41/PREREGISTRATION_v4.1.md)). Test set: 2,000 plasmids never used during development, 3,337 alignment-checked pairs.

- **Same-lineage identification (L5):** F1 0.86 (95% CI 0.79 to 0.92; precision 0.79, recall 0.94), against 0.58 for MOB-suite secondary clusters, 0.34 for pling and below 0.01 for mge-cluster. Difference from MOB-suite: +0.29 (95% CI 0.17 to 0.42).
- **Related-backbone grouping (L1):** F1 0.38, against 0.26 for MOB-suite primary clusters. The confidence interval was too wide to show non-inferiority (pre-registered criterion not met).
- **Stability:** 100% of codes unchanged in all 15 database-growth scenarios.
- **Reproducibility:** 100% of 500 plasmids re-typed from sequence received their published code.
- **Robustness:** 100% kept L1-L5 after 0.1% random substitutions (97% kept L6).
- **Speed (8 threads):** instant for a plasmid already in the database; 0.25 s per new plasmid related to the database and 1.08 s when divergent, against 1.25 s for MOB-suite on the same machine.
- **Replicon prediction:** 91.0% cross-validated accuracy, macro-F1 0.67, over 29 groups (8,404 training plasmids). Accuracy is high mainly for common groups.

All numbers are generated from the result files by [`build_facts.py`](build_facts.py) into [`docs/PLIN_FACTS.json`](docs/PLIN_FACTS.json).

### Database

Release **db-2026.10.05**: 127,517 unique plasmids, 93,440 L6 codes, 72,685 lineages (L5) and 27,333 backbone families (L1).

---

## Installation

### Desktop app (recommended: no Python or conda needed)

1. Go to the [latest release](https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification/releases/latest) and download the zip for your system:
   `pLIN-windows.zip` (Windows 10/11), `pLIN-macos.zip` (Apple silicon Macs) or `pLIN-linux.zip` (64-bit Linux from 2024 or later).
2. Unzip it and open **`README_FIRST.txt`** inside: it gives the exact steps for your system.
   In short:
   - **Windows:** Extract All, open the `pLIN` folder, double-click `pLIN.exe` (at "Windows protected your PC": More info, Run anyway).
   - **macOS:** drag `pLIN.app` to Applications; the first time, right-click it and choose Open, then Open.
   - **Linux:** `chmod +x pLIN/pLIN` and run `./pLIN/pLIN`.
3. pLIN opens in your web browser. Keep the code scheme on **v4.1** and click **Download the pLIN v4.1 database (1.6 GB)** (once).
4. Try the example: upload the 8 files in `sample_data/swiss_vim1_outbreak` (included in the zip), click **Run Analysis** and compare with `expected_pLIN_results.tsv`. The first analysis builds a search index once (12 to 16 GB, several minutes).
5. To stop pLIN, click **Quit pLIN** in the left sidebar.

Allow about 20 GB of free disk space. Full instructions and troubleshooting: [Desktop App User Manual](desktop_app/USER_MANUAL.md).

### From source (for developers and command-line use)

The steps below install pLIN with Python and conda.

### Prerequisites

- **Python 3.10 or higher** (Python 3.11+ recommended)
- **Git** (for cloning the repository)
- **Conda** (recommended) or **pip** with virtual environment
- **MMseqs2** (for typing new plasmids with pLIN v4.1): `conda install -c bioconda mmseqs2`

### Quick Start (All Platforms)

```bash
# 1. Clone the repository
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
cd pLIN-plasmid-classification

# 2. Run the install script (see platform-specific instructions below)
# macOS/Linux:
bash install_pLIN.sh

# Windows:
install_pLIN.bat

# 3. Download the pLIN v4.1 database (see below) into data/plin_v41/

# 4. Launch the GUI
conda activate pLIN_tools
streamlit run plin_app.py
```

---

### Database version & updates

pLIN separates two version numbers that change on different schedules:

| | Tracks | Where to check |
|---|---|---|
| **App version** | Code: the Streamlit app, classification logic, modules | `PLIN_APP_VERSION` in `plin_app.py`; shown in the app's Export tab |
| **Database version** (currently `db-2026.10.05`) | Content: the plasmid codes, protein-family catalogue and k-mer index | `data/plin_v41/DATABASE_VERSION.json` |

**Download the database.** The pLIN v4.1 database is about 1.6 GB, too large for the code repository. Download the files of release `db-2026.10.05` from the [Releases](../../releases) page or from its permanent archive on Zenodo ([doi:10.5281/zenodo.23126057](https://doi.org/10.5281/zenodo.23126057)) and place them in `data/plin_v41/`:

```
plin_v41_codes.tsv.gz        codes for every database plasmid (also usable on its own)
plin_v41_index.npz           codes, k-mer sketches and the backbone founder tree
family_members.faa.gz        protein-family catalogue (MMseqs2 search target)
family_exact_lookup.npz      protein-sequence hash -> family
plasmid_hashes.tsv.gz        whole-plasmid hash -> accession
DATABASE_VERSION.json        counts, levels, software versions, SHA-256 checksums
```

The desktop app and the Streamlit app can also download and check the files for you (a **Download the pLIN v4.1 database** button appears when the database is missing; files go to `~/.plin/plin_v41`). Check the files against the SHA-256 checksums in `DATABASE_VERSION.json`. On first use, the app builds an MMseqs2 search index of the protein catalogue (about 16 GB, a few minutes) in `~/.plin/v41_mmseqs_index`; set `PLIN_V41_INDEX` to put it elsewhere.

**Reproducibility note.** A plasmid already in the database always gets its published code. A new plasmid gets the code it would receive in the next release; levels that do not exist in the release are reported as provisional. Plasmids analysed together are coded consistently with each other. Codes become permanent for everyone once a plasmid is added to an official database release.

The earlier v3 codes (4-mer composition) remain available in the app as a legacy option and are listed next to the v4.1 codes in `plin_v41_codes.tsv.gz`.

---

### macOS Installation (Step-by-Step)

#### Step 1: Install Homebrew (if not installed)
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

#### Step 2: Install Python and Conda
```bash
# Option A: Install Miniconda (recommended)
brew install --cask miniconda
conda init zsh  # or: conda init bash
# Restart your terminal after this

# Option B: Install Python directly
brew install python@3.11
```

#### Step 3: Create Conda Environment
```bash
conda create -n pLIN_tools python=3.11 -y
conda activate pLIN_tools
```

#### Step 4: Install Python Dependencies
```bash
pip install -r requirements.txt
```

#### Step 5: Install MMseqs2 (required for pLIN v4.1) and optional tools
```bash
# MMseqs2: matches the proteins of new plasmids to the pLIN protein families (required)
conda install -c bioconda -c conda-forge "mmseqs2=18.8cc5c" -y

# AMRFinderPlus (AMR gene detection)
conda install -c bioconda -c conda-forge ncbi-amrfinderplus -y
amrfinder --update  # Download latest database

# MOBsuite (mobility typing)
pip install mob_suite

# Mash (MinHash ANI estimation)
conda install -c bioconda mash -y

# FastANI (true ANI computation)
conda install -c bioconda fastani -y

# minimap2 (SNP sub-typing)
conda install -c bioconda minimap2 -y

# MinCED (CRISPR spacer extraction)
conda install -c bioconda minced -y

# BLAST+ (CRISPR host inference)
conda install -c bioconda blast -y

# Prodigal (gene annotation)
conda install -c bioconda prodigal -y

# Ollama (AI chatbot: optional)
brew install ollama
ollama pull llama3.2
```

#### Step 6: Launch pLIN
```bash
conda activate pLIN_tools
streamlit run plin_app.py
```
The app will open automatically at `http://localhost:8501`.

---

### Linux Installation (Ubuntu/Debian: Step-by-Step)

#### Step 1: Install System Dependencies
```bash
sudo apt update
sudo apt install -y python3 python3-pip python3-venv git wget curl
```

#### Step 2: Install Miniconda
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda init bash
# Restart your terminal
```

#### Step 3: Create Conda Environment
```bash
conda create -n pLIN_tools python=3.11 -y
conda activate pLIN_tools
```

#### Step 4: Install Python Dependencies
```bash
pip install -r requirements.txt
```

#### Step 5: Install MMseqs2 (required for pLIN v4.1) and optional tools
```bash
# MMseqs2: matches the proteins of new plasmids to the pLIN protein families (required)
conda install -c bioconda -c conda-forge "mmseqs2=18.8cc5c" -y

# AMRFinderPlus
conda install -c bioconda -c conda-forge ncbi-amrfinderplus -y
amrfinder --update

# MOBsuite
pip install mob_suite

# Mash, FastANI, minimap2, MinCED, BLAST+, Prodigal
conda install -c bioconda mash fastani minimap2 minced blast prodigal -y

# Ollama (AI chatbot: optional)
curl -fsSL https://ollama.com/install.sh | sh
ollama pull llama3.2
```

#### Step 6: Launch pLIN
```bash
conda activate pLIN_tools
streamlit run plin_app.py
```

---

### Windows Installation (Step-by-Step)

#### Step 1: Install Miniconda
1. Download Miniconda from: https://docs.conda.io/en/latest/miniconda.html
2. Run the installer (`Miniconda3-latest-Windows-x86_64.exe`)
3. Check "Add Miniconda3 to my PATH" during installation
4. Open **Anaconda Prompt** from the Start Menu

#### Step 2: Create Conda Environment
```cmd
conda create -n pLIN_tools python=3.11 -y
conda activate pLIN_tools
```

#### Step 3: Install Python Dependencies
```cmd
pip install -r requirements.txt
```

#### Step 4: Install MMseqs2 (required for pLIN v4.1) and optional tools
```cmd
:: MMseqs2 (required): bioconda has no Windows build, so download the official Windows build
:: into the pLIN folder, where pLIN finds it (install_pLIN.bat does this for you)
powershell -Command "Invoke-WebRequest https://github.com/soedinglab/MMseqs2/releases/download/18-8cc5c/mmseqs-win64.zip -OutFile mmseqs-win64.zip; Expand-Archive mmseqs-win64.zip ."

:: AMRFinderPlus
conda install -c bioconda -c conda-forge ncbi-amrfinderplus -y
amrfinder --update

:: Mash, FastANI, minimap2 (via conda)
conda install -c bioconda mash fastani minimap2 minced blast prodigal -y

:: MOBsuite
pip install mob_suite

:: Ollama (download from https://ollama.com/download/windows)
:: After installing, run: ollama pull llama3.2
```

#### Step 5: Launch pLIN
```cmd
conda activate pLIN_tools
streamlit run plin_app.py
```

> **Note:** Some bioinformatics tools (AMRFinderPlus, MOBsuite) have limited Windows support. For full functionality, consider using **WSL2** (Windows Subsystem for Linux) and following the Linux instructions.

---

### Docker Installation (All Platforms)

```bash
# Build the Docker image
docker build -t plin .

# Run the container
docker run -p 8501:8501 -v $(pwd)/data:/app/data plin

# Access at http://localhost:8501
```

---

## Usage Guide

### Basic Workflow

1. **Launch the app:** `streamlit run plin_app.py`
2. **Upload FASTA files:** Drag and drop plasmid sequences (.fasta, .fa, .fna)
3. **Configure analysis:**
   - Select Inc group (or use auto-detect)
   - Enable/disable AMRFinderPlus, MOBsuite, Prodigal
   - Upload metadata CSV (optional: for temporal outbreak analysis)
   - Enable Mash ANI, FastANI, SNP sub-typing (optional)
4. **Click "Run pLIN Analysis"**
5. **Explore results** across 8 tabs

### Application Tabs

| Tab | Description |
|-----|-------------|
| **Overview** | pLIN system description, threshold table, sequence length warnings, metrics dashboard |
| **Results** | Interactive data table with search/filter, pLIN distribution, Inc group breakdown |
| **Cladogram** | Rectangular, circular, heatmap, and AMR-annotated cladograms |
| **AMR Analysis** | Gene prevalence, drug class pie charts, critical gene alerts, heatmaps |
| **Epidemiology** | Mobility prediction, outbreak detection, temporal clusters, Mash/FastANI, SNP sub-typing |
| **CRISPR Host** | CRISPR spacer extraction, host-plasmid heatmap, probability ranking |
| **Bacterial Buddy** | AI chatbot (Ollama LLM) for context-aware Q&A about your analysis |
| **Export** | Download TSV tables, PNG/PDF figures, ZIP bundle |

### Uploading Metadata

To enable temporal outbreak clustering and epidemiological analysis:

1. Prepare a CSV or TSV file with columns such as:
   - `plasmid_id` or `sample_id` (to match with FASTA files)
   - `collection_date` (any standard date format)
   - `location`, `hospital`, `ward` (optional)
2. Upload via the "Metadata CSV/TSV" uploader in the sidebar
3. Date columns are auto-detected and parsed
4. Metadata is merged into the integrated results table

### Optional Tool Integration

All external tools are **optional**: pLIN works without them but gains additional features when they are available:

| Tool | Feature Enabled | Install Command |
|------|----------------|-----------------|
| AMRFinderPlus | AMR/stress/virulence gene detection | `conda install -c bioconda ncbi-amrfinderplus` |
| MOBsuite | Relaxase family + MPF type classification | `pip install mob_suite` |
| Mash | Fast MinHash ANI estimation | `conda install -c bioconda mash` |
| FastANI | True average nucleotide identity | `conda install -c bioconda fastani` |
| minimap2 | SNP sub-typing within L6 clusters | `conda install -c bioconda minimap2` |
| MinCED | CRISPR spacer extraction | `conda install -c bioconda minced` |
| BLAST+ | CRISPR host inference | `conda install -c bioconda blast` |
| Prodigal | Gene/ORF annotation | `conda install -c bioconda prodigal` |
| Ollama | AI chatbot (Bacterial Buddy) | See platform-specific instructions |

---

## pLIN Classification System

The six levels and their thresholds are listed in the [Overview](#overview). Reading a code from the Swiss VIM-1 sample data:

```
169.178.183.208.209.2438
|   |   |   |   |   +-- L6: near-identical clone (k-mer similarity >= 0.95)
|   |   |   |   +------ L5: lineage (k-mer similarity >= 0.80)
|   |   |   +---------- L4: backbone variant (k-mer containment >= 0.80)
|   |   +-------------- L3: shared backbone (k-mer containment >= 0.50)
|   +------------------ L2: backbone group (protein-family containment >= 0.60)
+---------------------- L1: backbone family (protein-family containment >= 0.40)
```

Two plasmids sharing the first five numbers belong to the same lineage; sharing all six means they are near-identical. Numbers are identifiers, not distances: 209 and 210 are not more related than 209 and 900.

---

## Project Structure

```
pLIN-plasmid-classification/
├── plin_app.py                    # Main Streamlit GUI application
├── assign_pLIN.py                 # Batch pLIN assignment script
├── assign_pLIN_reference.py       # Reference database pLIN assignment
├── build_inc_centroids.py         # Train Inc group classifier
├── integrate_pLIN_AMR.py          # Merge pLIN + AMRFinderPlus results
├── install_pLIN.sh                # macOS/Linux install script
├── install_pLIN.bat               # Windows install script
├── launch_pLIN.command            # macOS double-click launcher
├── launch_pLIN.sh                 # Linux launcher
├── launch_pLIN.bat                # Windows launcher
├── requirements.txt               # Python dependencies
├── Dockerfile                     # Docker deployment
├── LICENSE                        # GPL-3.0 license
├── CITATION.cff                   # Citation metadata
├── plin_v41_typer.py              # pLIN v4.1 typing engine (used by the app and CLI)
├── plin_v41_nn.py, plin_v41.py    # v4.1 founder tree and nearest-relative assignment
├── plin_kmers.py                  # k-mer sketches (FracMinHash, k = 21)
├── build_facts.py                 # Generates docs/PLIN_FACTS.json from the result files
├── docs/PLIN_FACTS.json           # Every published number, with its source file
├── data/
│   ├── plin_v41/                  # v4.1 database (download from Releases)
│   ├── inc_classifier.npz         # Replicon (Inc/Rep) classifier (29 groups, 8,404 training plasmids)
│   └── inc_centroids.npz          # Inc group centroids
├── test_plasmids/
│   └── IncX/ (22 test FASTA files)
└── output/
    ├── pLIN_assignments.tsv              # v3 (legacy) codes, training plasmids
    ├── pLIN_reference_assignments.tsv     # v3 (legacy) codes, db-2026.10.02
    ├── reference_inc_classifications.tsv  # KNN Inc type classifications
    ├── integrated/
    │   ├── pLIN_AMR_integrated.tsv       # pLIN + AMRFinderPlus merged table
    │   └── pLIN_lineage_AMR_summary.tsv  # Lineage-level AMR summaries
    └── amrfinder/
        └── amrfinder_all_plasmids.tsv    # Raw AMRFinderPlus output (64,891 detections)
```

---

## Replicon (Inc/Rep) Groups Supported (29)

ColE, ColRNAI, IncA, IncAC2, IncC, IncF, IncFIB, IncFIBK, IncFIC, IncFII, IncHI1, IncHI2, IncI, IncI1, IncI2, IncLM, IncN, IncR, IncX1, IncX3, IncX4, repAci1, repAci_large, repEF_conj, repEF_res, repPae_large, repPae_small, repSA_large, repSA_small

---

## Troubleshooting

### Common Issues

| Issue | Solution |
|-------|----------|
| `ModuleNotFoundError` | Activate conda env: `conda activate pLIN_tools` |
| AMRFinderPlus not found | Install: `conda install -c bioconda ncbi-amrfinderplus` then `amrfinder --update` |
| Streamlit won't start | Check port: `streamlit run plin_app.py --server.port 8502` |
| Memory error (large dataset) | Reduce number of input files or increase system RAM |
| Mash/FastANI not detected | Install via conda and ensure PATH is set |
| Ollama connection error | Start Ollama: `ollama serve` then `ollama pull llama3.2` |

### Checking Tool Availability

The pLIN app auto-detects all optional tools at startup. Check the sidebar for tool status indicators.

---

## Citation

If you use pLIN in your research, you **must** cite:

> Xavier BB, Bari AK, Sinha B, Rossen JWA. pLIN: a permanent, multi-resolution nomenclature for bacterial plasmids, from backbone families to outbreak clones. Manuscript under review. https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification

Earlier version (preprint):

> Xavier BB, Bari AK, Sinha B, Rossen JWA. Development and validation of pLIN, a permanent lineage-numbering system for tracking antimicrobial resistance plasmids. Research Square (2026). https://doi.org/10.21203/rs.3.rs-10481391/v1

---

## License

Copyright (C) 2025-2026 Basil Britto Xavier. This project is licensed under the **GNU General Public License v3.0 or later** (see [LICENSE](LICENSE)). If you use pLIN in published work, please cite it as described in [CITATION.cff](CITATION.cff).
