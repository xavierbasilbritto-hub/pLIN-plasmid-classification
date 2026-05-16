# pLIN — Plasmid Lineage Identification Number

A hierarchical, reference-free classification system for bacterial plasmid genomes with integrated antimicrobial resistance (AMR) gene surveillance.

**Author:** Basil Xavier Britto | **License:** GPL-3.0 with mandatory citation clause | **Citation Required:** See [CITATION.cff](CITATION.cff)

---

## Overview

pLIN assigns each plasmid a **six-position hierarchical code** (e.g., `1.1.3.5.12.45`) based on tetranucleotide (4-mer) composition distances and single-linkage clustering at six biologically calibrated thresholds. The system spans from broad family-level (~85% ANI) to strain-level (~99.9% ANI) resolution.

### Key Features

| Category | Features |
|----------|----------|
| **Classification** | 6-level hierarchical pLIN codes (L1-L6), KNN Inc/Rep group detection (91.1% accuracy, 28 groups), Unknown/Novel flagging |
| **AMR Surveillance** | AMRFinderPlus integration (AMR + stress + virulence genes), critical gene alerts, drug class analysis |
| **Genomic Analysis** | Mash/MinHash ANI estimation, FastANI true ANI, minimap2 SNP sub-typing within L6 clusters |
| **Epidemiology** | Plasmid mobility prediction (MOBsuite + AMRFinderPlus), outbreak detection, temporal outbreak clustering (30-day window) |
| **Host Inference** | CRISPR spacer-based host prediction (MinCED + BLAST+), reference DB (72,959 plasmids) |
| **AI/ML** | Nucleotide Transformer LLM (optional), Bacterial Buddy AI chatbot (Ollama), adaptive per-Inc thresholds |
| **Visualization** | Interactive Streamlit GUI (8 tabs), cladograms, heatmaps, Plotly charts |
| **Deployment** | Cross-platform (macOS/Windows/Linux), Docker support, one-click launchers |

### Performance Metrics

- **Simpson's Diversity Index:** 0.985
- **Inc/Rep Detection Accuracy:** 91.1% (5-fold CV, 28 groups)
- **Training Dataset:** 8,077 plasmid sequences across 28 Inc/Rep groups (8,056 unique plasmids)
- **Reference Database:** 79,305 plasmids (8,056 training + 71,249 PLSDB/NCBI RefSeq)
- **Unique pLIN Codes:** 57,886 strain-level codes (across 79,305 plasmids)
- **Processing Time:** <30 minutes on a standard laptop

---

## Installation

### Prerequisites

- **Python 3.10 or higher** (Python 3.11+ recommended)
- **Git** (for cloning the repository)
- **Conda** (recommended) or **pip** with virtual environment

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

# 3. Launch the GUI
conda activate pLIN_tools
streamlit run plin_app.py
```

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
pip install streamlit numpy pandas scipy biopython scikit-learn matplotlib seaborn plotly python-pptx requests
```

#### Step 5: Install Optional Bioinformatics Tools
```bash
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

# Ollama (AI chatbot — optional)
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

### Linux Installation (Ubuntu/Debian — Step-by-Step)

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
pip install streamlit numpy pandas scipy biopython scikit-learn matplotlib seaborn plotly python-pptx requests
```

#### Step 5: Install Optional Bioinformatics Tools
```bash
# AMRFinderPlus
conda install -c bioconda -c conda-forge ncbi-amrfinderplus -y
amrfinder --update

# MOBsuite
pip install mob_suite

# Mash, FastANI, minimap2, MinCED, BLAST+, Prodigal
conda install -c bioconda mash fastani minimap2 minced blast prodigal -y

# Ollama (AI chatbot — optional)
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
pip install streamlit numpy pandas scipy biopython scikit-learn matplotlib seaborn plotly python-pptx requests
```

#### Step 4: Install Optional Bioinformatics Tools
```cmd
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
   - Upload metadata CSV (optional — for temporal outbreak analysis)
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

All external tools are **optional** — pLIN works without them but gains additional features when they are available:

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

### Hierarchical Levels

| Level | Bin | Cosine Distance (d) | ANI Equivalent | Biological Meaning |
|-------|-----|---------------------|----------------|-------------------|
| L1 | A | d <= 0.150 | ~85% | Broad plasmid family |
| L2 | B | d <= 0.100 | ~90% | Subfamily |
| L3 | C | d <= 0.050 | ~95% | Cluster (species-level) |
| L4 | D | d <= 0.020 | ~98% | Subcluster |
| L5 | E | d <= 0.010 | ~99% | Clone complex |
| L6 | F | d <= 0.001 | ~99.9% | Strain / Outbreak |

### Example pLIN Code
```
1.1.3.5.12.45
| | | | |  +-- L6: Strain-level cluster (d <= 0.001)
| | | | +---- L5: Clone complex (d <= 0.010)
| | | +------ L4: Subcluster (d <= 0.020)
| | +-------- L3: Cluster (d <= 0.050)
| +---------- L2: Subfamily (d <= 0.100)
+------------ L1: Family (d <= 0.150)
```

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
├── data/
│   ├── inc_classifier.npz         # Trained KNN classifier (28 groups, 8,077 samples)
│   └── inc_centroids.npz          # Inc group centroids
├── test_plasmids/
│   └── IncX/ (22 test FASTA files)
└── output/
    ├── pLIN_assignments.tsv              # 8,056 training plasmid assignments
    ├── pLIN_reference_assignments.tsv     # 79,305 full reference database assignments
    ├── reference_inc_classifications.tsv  # KNN Inc type classifications
    ├── integrated/
    │   ├── pLIN_AMR_integrated.tsv       # pLIN + AMRFinderPlus merged table
    │   └── pLIN_lineage_AMR_summary.tsv  # Lineage-level AMR summaries
    └── amrfinder/
        └── amrfinder_all_plasmids.tsv    # Raw AMRFinderPlus output (64,891 detections)
```

---

## Inc Groups Supported (20)

ColE, ColRNAI, IncA, IncAC2, IncC, IncF, IncFIB, IncFIBK, IncFIC, IncFII, IncHI1, IncHI2, IncI, IncI1, IncI2, IncN, IncR, IncX1, IncX3, IncX4

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

> Xavier, B. (2025). pLIN: A Plasmid Lineage Identification Number System for Hierarchical, Permanent Classification of Bacterial Plasmids Integrated with Antimicrobial Resistance Gene Surveillance. https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification

---

## License

This project is licensed under **GPL-3.0** with a mandatory citation clause. See [LICENSE](LICENSE) and [CITATION.cff](CITATION.cff) for details.
