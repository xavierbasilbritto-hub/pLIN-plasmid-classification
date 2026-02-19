# pLIN Installation Guide

**pLIN v2.1.0 — Plasmid Life Identification Number System**

Three installation methods are available. Choose based on your preference:

| Method | Best for | Time | Prerequisites |
|--------|----------|------|---------------|
| **A. Docker (Recommended)** | Reviewers, quick testing | 5-10 min | Docker Desktop only |
| **B. Automated script** | All platforms | 5-15 min | Python 3.10+ |
| **C. Manual install** | Full control | 10-20 min | Python 3.10+, conda |

---

## Method A: Docker (Recommended for Reviewers)

Docker provides a fully self-contained environment with all tools pre-installed. No dependency conflicts possible.

### Prerequisites

Install Docker Desktop for your operating system:
- **macOS**: https://docs.docker.com/desktop/install/mac-install/
- **Windows**: https://docs.docker.com/desktop/install/windows-install/
- **Linux**: https://docs.docker.com/desktop/install/linux-install/

### Step 1: Clone the repository

```bash
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
cd pLIN-plasmid-classification
```

### Step 2: Build and run

```bash
docker compose up --build
```

Or equivalently:

```bash
docker build -t plin .
docker run -p 8501:8501 -v ./data:/app/data -v ./output:/app/output plin
```

### Step 3: Open in browser

Navigate to: **http://localhost:8501**

### Step 4: Stop

Press `Ctrl+C` in the terminal, or run:

```bash
docker compose down
```

### Troubleshooting (Docker)

| Issue | Solution |
|-------|----------|
| Port 8501 in use | Change port: `docker run -p 9501:8501 plin` then open http://localhost:9501 |
| Build fails on ARM/Apple Silicon | The image auto-detects architecture; if issues persist, try: `docker build --platform linux/amd64 -t plin .` |
| Out of memory | Increase Docker Desktop memory limit to 4 GB (Settings > Resources) |

---

## Method B: Automated Setup Script (All Platforms)

A single Python script handles installation and launch on any OS.

### macOS

```bash
# 1. Install Python 3.10+ (if not already installed)
brew install python@3.11

# 2. Clone the repository
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
cd pLIN-plasmid-classification

# 3. Run the setup script
python3 setup_pLIN.py
```

### Linux (Ubuntu/Debian)

```bash
# 1. Install Python 3.10+ (if not already installed)
sudo apt update
sudo apt install python3 python3-pip python3-venv git

# 2. Clone the repository
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
cd pLIN-plasmid-classification

# 3. Run the setup script
python3 setup_pLIN.py
```

### Linux (Fedora/RHEL)

```bash
# 1. Install Python 3.10+
sudo dnf install python3 python3-pip git

# 2. Clone and run
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
cd pLIN-plasmid-classification
python3 setup_pLIN.py
```

### Windows

```powershell
# 1. Install Python 3.10+ from https://www.python.org/downloads/
#    IMPORTANT: Check "Add Python to PATH" during installation

# 2. Open Command Prompt or PowerShell, then:
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
cd pLIN-plasmid-classification

# 3. Run the setup script
python setup_pLIN.py
```

### Setup script options

```bash
python3 setup_pLIN.py                # Full install + launch
python3 setup_pLIN.py --install      # Install only (no launch)
python3 setup_pLIN.py --launch       # Launch only (skip install)
python3 setup_pLIN.py --check        # Check environment status
python3 setup_pLIN.py --no-biotools  # Skip optional bioinformatics tools
python3 setup_pLIN.py --docker       # Build and run via Docker
```

---

## Method C: Manual Installation

### Step 1: Install Python 3.10+

| OS | Command |
|----|---------|
| macOS | `brew install python@3.11` |
| Ubuntu/Debian | `sudo apt install python3 python3-pip python3-venv` |
| Fedora | `sudo dnf install python3 python3-pip` |
| Windows | Download from https://www.python.org/downloads/ |

### Step 2: Clone the repository

```bash
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
cd pLIN-plasmid-classification
```

### Step 3: Create virtual environment

```bash
# macOS / Linux
python3 -m venv .venv
source .venv/bin/activate

# Windows
python -m venv .venv
.venv\Scripts\activate
```

### Step 4: Install Python dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

### Step 5: Launch pLIN

```bash
streamlit run plin_app.py
```

The app opens at **http://localhost:8501**.

### Step 6 (Optional): Install bioinformatics tools

These are optional. Without them, core pLIN classification and pLIN assignment work. With them, additional features are available (AMR detection, CRISPR analysis, ANI validation).

```bash
# Install Miniconda first if not already installed:
# https://docs.conda.io/en/latest/miniconda.html

# Configure channels
conda config --add channels bioconda
conda config --add channels conda-forge

# Install tools
conda install -c bioconda -c conda-forge \
    ncbi-amrfinderplus mash fastani minimap2 minced blast prodigal

# Update AMRFinderPlus database
amrfinder --update
```

---

## Quick Test After Installation

1. Open **http://localhost:8501** in your browser
2. The pLIN GUI should display with the sidebar and main analysis panel
3. Upload a FASTA file (any plasmid sequence) to test pLIN assignment
4. The tool will:
   - Compute 4-mer frequency vectors
   - Classify the Inc group via KNN
   - Assign a hierarchical pLIN code via nearest-neighbour lookup
   - Display results with nearest-neighbour distance and confidence

### Sample test data

Test FASTA files can be obtained from NCBI:

```bash
# Download a sample IncN plasmid (pKPC-2)
curl -o test_plasmid.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NZ_CP015070.1&rettype=fasta"
```

Upload `test_plasmid.fasta` to the pLIN app to verify the installation.

---

## Feature Availability by Installation Method

| Feature | Docker | Script (core) | Script (+ conda) | Manual (core) | Manual (+ conda) |
|---------|:------:|:-------------:|:-----------------:|:-------------:|:-----------------:|
| pLIN classification | Yes | Yes | Yes | Yes | Yes |
| Inc group prediction | Yes | Yes | Yes | Yes | Yes |
| Single-plasmid query mode | Yes | Yes | Yes | Yes | Yes |
| Multi-plasmid clustering | Yes | Yes | Yes | Yes | Yes |
| AMR gene detection | Yes | No | Yes | No | Yes |
| CRISPR spacer analysis | Yes | No | Yes | No | Yes |
| Mash/FastANI validation | Yes | No | Yes | No | Yes |
| SNP sub-typing | Yes | No | Yes | No | Yes |
| Outbreak detection | Yes | Yes | Yes | Yes | Yes |

---

## System Requirements

| Requirement | Minimum | Recommended |
|-------------|---------|-------------|
| RAM | 2 GB | 4 GB |
| Disk space | 50 MB (core) | 500 MB (with tools) |
| Python | 3.10 | 3.11 |
| Docker | 20.10+ | Latest |
| Browser | Any modern browser | Chrome/Firefox |

---

## Contact

For issues or questions:
- GitHub Issues: https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification/issues
- Correspondence: [corresponding author email]
