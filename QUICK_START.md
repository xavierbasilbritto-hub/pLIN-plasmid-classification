# pLIN Quick Start Guide

**For reviewers with no bioinformatics experience.**
**Time needed: 20-30 minutes (mostly waiting for downloads).**
**Works on: macOS, Linux, and Windows.**

---

## What is pLIN?

pLIN is a web application that runs locally on your computer. You open it in your web browser (Chrome, Safari, Firefox, or Edge) — just like any website, except nothing is uploaded to the internet.

You give it plasmid DNA sequence files, and it:
- Classifies them into Inc/Rep groups (28 supported)
- Assigns a 6-level hierarchical code (like `3.5.12.45.201.3050`)
- Detects antimicrobial resistance (AMR), virulence, and stress genes
- Filters chromosomal contigs from plasmid sequences automatically
- Detects IS elements (insertion sequences) and mobile genetic elements
- Shows interactive trees, charts, and gene maps

---

## Before You Start

| Requirement | macOS | Linux | Windows |
|-------------|-------|-------|---------|
| **OS version** | macOS 13 (Ventura)+ | Ubuntu 20.04+, Debian 11+, Fedora 36+, CentOS 8+, Arch, RHEL 8+ | Windows 10 (1903+) or 11 |
| **Disk space** | 5 GB | 5 GB | 5 GB |
| **Admin access** | No | sudo access | No |
| **Internet** | Setup only | Setup only | Setup only |
| **Time** | 20-30 min | 20-30 min | 25-40 min |

---

## Step 1: Open a Terminal

### macOS
1. Press **Cmd + Space** (opens Spotlight search)
2. Type **Terminal**
3. Press **Enter**

> **Tip:** Paste in Terminal with **Cmd + V**.

### Linux
- **Ubuntu/Debian:** Press **Ctrl + Alt + T**
- **Fedora/RHEL:** Click "Activities" (top-left), type **Terminal**, click it
- **Any Linux:** Right-click the desktop and choose **"Open Terminal"**

> **Tip:** Paste in Terminal with **Ctrl + Shift + V**.

### Windows
- Open **Start Menu**, search for **Anaconda Prompt** (installed in Step 2)
- If you don't have Anaconda Prompt yet, use **Command Prompt** (press **Win + R**, type `cmd`, press Enter)

> **Tip:** Paste in Command Prompt by right-clicking.

---

## Step 2: Install Prerequisites

### macOS

**Install Homebrew** (if you don't have it):
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

**Install Miniconda:**
```bash
brew install --cask miniconda
conda init zsh
```
Close and reopen Terminal after this.

### Linux

**Install basic tools:**

Ubuntu / Debian:
```bash
sudo apt update && sudo apt install -y git curl wget unzip
```

Fedora / RHEL:
```bash
sudo dnf install -y git curl wget unzip
```

Arch Linux:
```bash
sudo pacman -Sy --noconfirm git curl wget unzip
```

**Install Miniconda:**
```bash
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -p ~/miniconda3
~/miniconda3/bin/conda init bash
```
Close and reopen Terminal after this.

### Windows

1. Download **Miniconda** from: https://docs.conda.io/en/latest/miniconda.html
2. Run the installer — **check "Add to PATH"** when asked
3. Click **Install**, then **Finish**
4. Open **Anaconda Prompt** from the Start Menu

---

## Step 3: Get the pLIN Package

If you received a ZIP file:
1. Unzip it (double-click on macOS/Windows, or `unzip` on Linux)
2. Navigate into the folder:

**macOS / Linux:**
```bash
cd ~/Downloads/PLASMID_TOOL
```

**Windows:**
```cmd
cd %USERPROFILE%\Downloads\PLASMID_TOOL
```

If cloning from GitHub:
```bash
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
cd pLIN-plasmid-classification
```

Verify you're in the right place — you should see `plin_app.py`:

```bash
ls plin_app.py
```

---

## Step 4: Run the Installer

### macOS / Linux
```bash
python3 setup_pLIN.py --install
```

### Windows
```cmd
python setup_pLIN.py --install
```

**What to expect** (15-30 minutes):
- `[OK]` messages for each step
- Conda environment creation + bioinformatics tool installation
- IS element reference database download (25 sequences from NCBI)
- BLAST database verification
- Some `[WARN]` messages for optional tools — this is normal

The installer will:
1. Check your Python version
2. Install system dependencies (Linux only)
3. Create a conda environment with bioinformatics tools
4. Install Python packages
5. Create output directories
6. Download IS reference sequences and build BLAST database
7. Verify BLAST availability

---

## Step 5: Launch pLIN

### macOS / Linux
```bash
python3 setup_pLIN.py --launch
```

### Windows
```cmd
python setup_pLIN.py --launch
```

Your browser will open automatically at **http://localhost:8501**.

---

## Step 6: Test with Sample Data

1. In the pLIN web interface, click **"Browse files"** or drag-and-drop
2. Navigate to `test_plasmids/IncX/` and upload 3 FASTA files
3. Click **"Run Analysis"**
4. Check the **Results** tab — you should see:
   - Inc/Rep group classification (e.g., IncX1, IncX3, IncX4)
   - 6-level pLIN codes (e.g., `3.5.12.45.201.3050`)
   - Confidence scores
5. Check the **Cladogram** tab for the interactive phylogenetic tree
6. Try the **Export** button to download a TSV file

---

## Step 7: Stop pLIN

Press **Ctrl + C** in the terminal window.

---

## Restarting pLIN Later

```bash
cd /path/to/PLASMID_TOOL
python3 setup_pLIN.py --launch
```

Or directly:
```bash
conda activate pLIN_tools
streamlit run plin_app.py
```

---

## Understanding pLIN Features

### Contig Classification (Plasmid vs Chromosome)

When you upload multi-contig assemblies or mixed FASTA files, pLIN automatically distinguishes plasmid sequences from chromosomal DNA before assigning pLIN codes. This prevents chromosomal contigs from receiving incorrect plasmid classifications.

**How it works:**

pLIN uses a multi-signal scoring system that combines four lines of evidence:

| Signal | What it checks | Score impact |
|--------|---------------|--------------|
| **Sequence length** | Sequences >500 kb are likely chromosomal | -50 to +30 |
| **KNN distance** | 4-mer (tetranucleotide) distance to 8,077 known plasmids | -25 to +30 |
| **Header keywords** | "plasmid", "chromosome", "genome" in FASTA headers | -20 to +15 |
| **Inc group confidence** | Strong Inc/Rep group match suggests plasmid | -5 to +15 |

**Three classification outcomes:**

1. **Plasmid** (high confidence, score >= 10) — receives a pLIN code and all downstream analyses (AMR, mobility, etc.)
2. **Incomplete plasmid** (ambiguous, -10 < score < 10) — gets AMR/mobility analysis but no pLIN code (prevents unreliable lineage assignments)
3. **Chromosome** (score <= -10) — excluded from all plasmid-specific analyses

**Multi-contig assembly handling:**

When multiple contigs from the same source file share the same Inc type, pLIN merges them into a single record (joined with 100-N spacer) before pLIN assignment. Contigs with different Inc types are kept separate as genuinely different plasmids.

**No external tools required** — contig classification uses the same KNN classifier (256-dimensional 4-mer frequency vectors + cosine distance) that powers Inc/Rep group detection. It runs automatically when "Auto-detect plasmid contigs" is enabled in the sidebar.

**Controlling the feature:**
- Enabled by default (recommended for mixed assemblies)
- Uncheck "Auto-detect plasmid contigs" in the sidebar to force pLIN assignment on all sequences (only if you are certain all inputs are plasmids)

---

### IS Element / Mobile Genetic Element (MGE) Detection

pLIN can detect insertion sequences (IS elements) in plasmid genomes using BLAST against a curated reference database of 25 IS families from ISfinder/NCBI.

**What IS elements are:**

Insertion sequences are short (800-2,500 bp) transposable elements that carry only their own transposase gene. They are the simplest autonomous transposons and are clinically important because:

- They mobilize antimicrobial resistance genes between plasmids and chromosomes
- Paired IS elements form **composite transposons** that capture and spread AMR genes
- Their distribution patterns reveal plasmid evolutionary relationships

**IS families detected (25):**

| Category | IS Families |
|----------|-------------|
| Gram-negative (18) | IS26, ISEcp1, IS1, IS903, IS6100, ISKpn26, IS5, IS3, IS4321, IS15, IS10, IS2, IS4, IS30, IS66, IS110, ISPa, ISAba |
| Gram-positive (7) | IS256, IS257/IS431, IS16, ISEnfa, IS1216, IS1251, Tn916 |

**How detection works:**

1. BLAST alignment against curated IS reference database
2. Filtering: >=80% identity, >=50% coverage of reference IS, e-value <=1e-10
3. Composite transposon detection: paired IS elements of the same family flanking AMR genes (identified from AMRFinderPlus results)

**Setup requirements:**

- **BLAST+** — installed automatically by `setup_pLIN.py`
- **IS reference database** — downloaded automatically during setup from NCBI
- If the database was not set up, run `python3 setup_pLIN.py --install` or `python3 detect_IS_elements.py`

**Where to find results:**

- In the pLIN GUI: **MGE Detection** tab (requires Prodigal gene predictions)
- Output directory: `output/mge_detection/`
  - `is_element_hits.tsv` — all IS hits per plasmid
  - `is_family_counts.tsv` — IS family distribution across groups
  - `composite_transposons.tsv` — IS pairs flanking AMR genes

---

## Troubleshooting

### All Platforms

**"ModuleNotFoundError: No module named 'streamlit'"**
```bash
conda activate pLIN_tools
pip install -r requirements.txt
```

**"Port 8501 is already in use"**
```bash
streamlit run plin_app.py --server.port 8502
```
Then open http://localhost:8502.

**Browser shows a blank page**
1. Wait 10 seconds, then refresh
2. Try http://127.0.0.1:8501 instead
3. Check terminal for error messages

**AMR Analysis tab says "tool not available"**
```bash
conda activate pLIN_tools
conda install -c bioconda -c conda-forge ncbi-amrfinderplus -y
amrfinder --update
```

**IS element detection shows no results**
```bash
python3 setup_pLIN.py --install
```
Or manually:
```bash
python3 detect_IS_elements.py
```

### macOS

**"command not found: conda"**
```bash
source ~/miniconda3/bin/activate
conda init zsh
```
Close and reopen Terminal.

### Linux

**"Permission denied"**
```bash
chmod +x setup_pLIN.py
```

**"E: Unable to locate package"**
```bash
sudo apt update
```

### Windows

**"conda is not recognized"**

Use **Anaconda Prompt** (not regular Command Prompt). Or reinstall Miniconda and check the "Add to PATH" box.

**Full bioinformatics tool support**

For complete BLAST/AMRFinder support, use WSL2:
```powershell
wsl --install
```
Then run pLIN inside the WSL Ubuntu terminal.

---

## Docker Alternative (All Platforms)

If you have Docker installed, you can skip the manual setup entirely:

### macOS / Linux
```bash
docker compose up --build
```

### Windows
Install Docker Desktop with WSL2 backend, then:
```cmd
docker compose up --build
```

Open **http://localhost:8501** in your browser. Stop with `docker compose down`.

---

## Check Your Environment

At any time, verify your installation:
```bash
python3 setup_pLIN.py --check
```

This shows the status of Python, packages, bioinformatics tools, conda environment, IS database, and BLAST.

---

## What Each File Does

| File | Purpose |
|------|---------|
| `plin_app.py` | Main application (do not edit) |
| `setup_pLIN.py` | Cross-platform installer and launcher |
| `requirements.txt` | Python package dependencies |
| `data/inc_classifier.npz` | KNN classifier (8,077 plasmids, 28 groups) |
| `data/inc_centroids.npz` | Group centroids for distance calculation |
| `output/pLIN_assignments.tsv` | Training database pLIN assignments |
| `output/pLIN_reference_assignments.tsv` | Full reference database (79,305 plasmids) |
| `output/mge_detection/` | IS element reference database and detection results |
| `detect_IS_elements.py` | Standalone IS element detection script |
| `build_IS_database.py` | Curated IS database builder |
| `validate_contig_classifier.py` | Contig classification validation suite |
| `test_plasmids/IncX/*.fasta` | Sample plasmid sequences for testing |
| `Dockerfile` | Docker-based installation |

---

## Platform-Specific Guides

For detailed, step-by-step instructions tailored to your operating system:

- [macOS Guide](QUICK_START_macOS.md) — Ventura, Sonoma, Sequoia (Intel & Apple Silicon)
- [Linux Guide](QUICK_START_Linux.md) — Ubuntu, Debian, Fedora, CentOS, Arch, RHEL, openSUSE
- [Windows Guide](QUICK_START_Windows.md) — Windows 10/11 with WSL2 option
