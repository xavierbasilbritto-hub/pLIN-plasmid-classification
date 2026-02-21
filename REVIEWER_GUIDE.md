# pLIN Reviewer Validation Guide

**Step-by-step instructions to install, run, and validate the pLIN tool.**
**Every step has been tested. Follow them exactly in order.**

---

## What Is pLIN?

pLIN (Plasmid Lineage Identification Number) is a web application that classifies bacterial plasmids. You upload plasmid DNA sequence files (FASTA format), and it:

- Assigns a hierarchical 6-level classification code (like `1.3.5.12.45.201`)
- Detects the incompatibility (Inc) group automatically
- Finds antimicrobial resistance (AMR) genes
- Detects outbreak-related plasmid clusters
- Produces interactive visualizations (trees, heatmaps, charts)

The tool runs as a local web app in your browser. Nothing is uploaded to the internet.

---

## Pick Your Operating System

| I am using... | Go to... |
|---|---|
| macOS (any Mac) | [Section A](#a-macos-setup) |
| Linux (Ubuntu, Debian, Fedora, etc.) | [Section B](#b-linux-setup) |
| Windows 10 or 11 | [Section C](#c-windows-setup) |
| Any OS with Docker installed | [Section D](#d-docker-setup-any-os) |

> **Recommendation for reviewers:** Docker (Section D) is the fastest and most reliable option. It works identically on all operating systems and requires no manual dependency installation.

---

## A. macOS Setup

### Prerequisites

- macOS 13 (Ventura) or newer
- At least 5 GB free disk space
- Internet connection (for first-time setup only)

### Step A1. Open Terminal

1. Press **Cmd + Space** to open Spotlight
2. Type **Terminal** and press Enter
3. A black/white window opens — this is your terminal

### Step A2. Install Homebrew (if you do not have it)

Paste this command and press Enter:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

- It will ask for your Mac password — type it (nothing shows while typing, that is normal)
- Wait until it says "Installation successful"
- **If you already have Homebrew:** skip this step

**Check it worked:**

```bash
brew --version
```

You should see something like `Homebrew 4.x.x`. If you see "command not found", follow the instructions Homebrew printed at the end of installation.

### Step A3. Install Miniconda

```bash
brew install --cask miniconda
```

Then initialize it for your terminal:

```bash
conda init zsh
```

**Close the terminal window and open a new one** (this is required).

**Check it worked:**

```bash
conda --version
```

You should see `conda 24.x.x` or similar.

### Step A4. Download pLIN

```bash
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
```

```bash
cd pLIN-plasmid-classification
```

**Check it worked:**

```bash
ls plin_app.py
```

You should see `plin_app.py` listed.

### Step A5. Run the Installer

```bash
chmod +x pLIN_macOS.sh
./pLIN_macOS.sh --install
```

This will:
- Create a conda environment called `pLIN_tools`
- Install Python packages (numpy, pandas, streamlit, etc.)
- Install bioinformatics tools (AMRFinderPlus, minimap2, BLAST+, etc.)

**This takes 15-30 minutes on first run.** You will see `[OK]` messages as each tool installs. Some tools may show `[WARN] Could not install` — that is fine, they are optional.

When finished, you will see:

```
  [OK]   Installation complete!
```

### Step A6. Launch pLIN

```bash
./pLIN_macOS.sh --launch
```

Your browser will automatically open to `http://localhost:8501`.

You should see:
- Title: **pLIN Classifier**
- A file upload area
- A sidebar with options

> **To stop the app:** go back to the terminal and press **Ctrl+C**.

**Now go to [Section E. Test the Tool](#e-test-the-tool).**

---

## B. Linux Setup

### Prerequisites

- Ubuntu 20.04+, Debian 11+, Fedora 36+, CentOS 8+, Arch Linux, or RHEL 8+
- At least 5 GB free disk space
- `sudo` access (for installing system packages)
- Internet connection

### Step B1. Open a Terminal

- Ubuntu/Debian: Press **Ctrl+Alt+T**
- Fedora/RHEL: Open "Activities" and search for "Terminal"
- Or right-click the desktop and select "Open Terminal"

### Step B2. Install Git (if you do not have it)

**Ubuntu/Debian:**

```bash
sudo apt update && sudo apt install -y git curl wget
```

**Fedora/RHEL:**

```bash
sudo dnf install -y git curl wget
```

**Arch:**

```bash
sudo pacman -Sy --noconfirm git curl wget
```

**Check it worked:**

```bash
git --version
```

You should see `git version 2.x.x`.

### Step B3. Download pLIN

```bash
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
```

```bash
cd pLIN-plasmid-classification
```

**Check it worked:**

```bash
ls plin_app.py
```

You should see `plin_app.py` listed.

### Step B4. Run the Installer

```bash
chmod +x pLIN_Linux.sh
./pLIN_Linux.sh --install
```

This will:
- Install system dependencies (build tools, Java for MinCED)
- Download and install Miniconda (if conda is not already installed)
- Create a conda environment called `pLIN_tools`
- Install all Python packages
- Install bioinformatics tools

**This takes 15-30 minutes on first run.** You may be asked for your sudo password.

When finished, you will see:

```
  [OK]   Installation complete!
```

### Step B5. Launch pLIN

```bash
./pLIN_Linux.sh --launch
```

Your browser will open to `http://localhost:8501`.

If the browser does not open automatically, manually open this URL in your browser:

```
http://localhost:8501
```

> **To stop the app:** go back to the terminal and press **Ctrl+C**.

**Now go to [Section E. Test the Tool](#e-test-the-tool).**

---

## C. Windows Setup

### Prerequisites

- Windows 10 (version 1903 or newer) or Windows 11
- At least 5 GB free disk space
- Internet connection

### Step C1. Install Git for Windows

1. Go to: https://git-scm.com/download/win
2. Download the installer and run it
3. Click "Next" on every screen (default settings are fine)
4. Click "Install"

### Step C2. Install Miniconda

1. Go to: https://docs.conda.io/en/latest/miniconda.html
2. Download **Miniconda3 Windows 64-bit**
3. Run the installer
4. **IMPORTANT:** On the "Advanced Options" screen, check the box that says **"Add Miniconda3 to my PATH environment variable"**
5. Click "Install"

### Step C3. Open Anaconda Prompt

1. Click the **Start Menu** (Windows icon, bottom-left)
2. Type **Anaconda Prompt**
3. Click **Anaconda Prompt** (not PowerShell, not cmd.exe)

> You must use Anaconda Prompt for all remaining steps.

### Step C4. Download pLIN

Type these commands in Anaconda Prompt:

```cmd
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
```

```cmd
cd pLIN-plasmid-classification
```

**Check it worked:**

```cmd
dir plin_app.py
```

You should see `plin_app.py` in the listing.

### Step C5. Run the Installer

```cmd
pLIN_Windows.bat --install
```

This will:
- Detect your Python version
- Create a conda environment called `pLIN_tools`
- Install all Python packages
- Attempt to install bioinformatics tools via conda

**This takes 15-30 minutes on first run.**

> **Note:** Some bioinformatics tools (AMRFinderPlus, Prodigal) have limited Windows support. For full functionality, consider using WSL2 (see step C7 below) or Docker (Section D).

When finished, you will see:

```
  Installation Complete!
```

### Step C6. Launch pLIN

```cmd
pLIN_Windows.bat --launch
```

Your browser will open to `http://localhost:8501`.

> **To stop the app:** go back to Anaconda Prompt and press **Ctrl+C**.

### Step C7. (Optional) Full Tool Support via WSL2

For all bioinformatics tools to work on Windows, install Windows Subsystem for Linux:

1. Open **PowerShell as Administrator** (right-click Start Menu > Terminal (Admin))
2. Run:
   ```powershell
   wsl --install
   ```
3. Restart your computer
4. Open the Ubuntu terminal from the Start Menu
5. Follow the **Linux Setup (Section B)** instructions inside WSL

**Now go to [Section E. Test the Tool](#e-test-the-tool).**

---

## D. Docker Setup (Any OS)

This is the easiest method. Docker runs pLIN inside a container with all tools pre-installed. No manual dependency installation required.

### Step D1. Install Docker

| Operating System | How to Install Docker |
|---|---|
| **macOS** | Download Docker Desktop from https://www.docker.com/products/docker-desktop/ — open the `.dmg` file and drag Docker to Applications. Launch Docker from Applications. |
| **Linux (Ubuntu/Debian)** | Run: `sudo apt update && sudo apt install -y docker.io docker-compose-plugin` then `sudo usermod -aG docker $USER` — then **log out and log back in**. |
| **Linux (Fedora/RHEL)** | Run: `sudo dnf install -y docker docker-compose-plugin` then `sudo systemctl start docker` then `sudo usermod -aG docker $USER` — then **log out and log back in**. |
| **Windows** | Download Docker Desktop from https://www.docker.com/products/docker-desktop/ — run the installer. When asked, enable WSL2 backend. Restart your computer. |

**Check it worked:**

```bash
docker --version
```

You should see `Docker version 24.x.x` or newer.

### Step D2. Download pLIN

```bash
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
```

```bash
cd pLIN-plasmid-classification
```

### Step D3. Build and Run pLIN

**Option 1: Docker Compose (recommended)**

```bash
docker compose up --build
```

**Option 2: Manual Docker commands**

```bash
docker build --build-arg INSTALL_BIOTOOLS=true -t plin-full .
```

```bash
docker run -p 8501:8501 -v ./data:/app/data -v ./output:/app/output plin-full
```

**The first build takes 10-20 minutes** (downloads Python packages + bioinformatics tools).

### Step D4. Open pLIN

Open your web browser and go to:

```
http://localhost:8501
```

You should see the pLIN Classifier interface.

### Step D5. Stop pLIN

**Docker Compose:**

```bash
docker compose down
```

**Manual Docker:**

```bash
docker stop plin-app
```

### Step D6. Run pLIN Again (After First Build)

**Docker Compose:**

```bash
docker compose up
```

**Manual Docker:**

```bash
docker run -p 8501:8501 -v ./data:/app/data -v ./output:/app/output plin-full
```

**Now go to [Section E. Test the Tool](#e-test-the-tool).**

---

## E. Test the Tool

This section walks you through a complete test using sample data included in the repository. Follow every step exactly.

### Step E1. Verify the App Is Running

Open your browser to:

```
http://localhost:8501
```

**You should see:**
- Page title: **pLIN Classifier**
- A file upload area that says "Upload one or more plasmid FASTA files"
- A sidebar on the left with configuration options

**If you see a blank page or error:**
- Go back to your terminal and check for error messages
- Make sure port 8501 is not used by another application
- Try: `http://127.0.0.1:8501` instead

### Step E2. Upload Test Plasmid Files

The repository includes 22 test plasmid sequences in the `test_plasmids/IncX/` folder.

1. In the pLIN web interface, find the **"Upload FASTA files"** area
2. Click **"Browse files"**
3. Navigate to the `pLIN-plasmid-classification/test_plasmids/IncX/` folder
4. Select **5 files** (for a quick test):
   - `IncX3_JN247852.fasta`
   - `SP12_P2.fasta`
   - `SP30_P3.fasta`
   - `SP42_P2.fasta`
   - `SP50_P4.fasta`
5. Click **Open**

**You should see:** "5 files uploaded" (or similar message) with the file names listed.

### Step E3. Configure the Analysis

In the **sidebar** (left side):

1. **Inc Group:** Leave as **"Auto-detect"** (default)
2. **Linkage Method:** Leave as **"Single"** (default)
3. **Adaptive Thresholds:** Leave **checked** (default)
4. If you see tool toggles (AMRFinderPlus, Prodigal, etc.):
   - Leave all **checked** if available
   - Tools that are not installed will be grayed out — this is normal

### Step E4. Run the Analysis

1. Click the **"Run Analysis"** button (blue button)
2. Wait for the analysis to complete
   - You will see a progress bar and status messages
   - This takes 10-30 seconds for 5 files

**When complete, you should see results appear in the main panel.**

### Step E5. Check the Results Tab

Click the **"Results"** tab.

**You should see a table with these columns:**
- `plasmid_id` — name of each plasmid
- `inc_type` — the detected incompatibility group (should show IncX variants)
- `pLIN` — the 6-level classification code (format: `X.X.X.X.X.X`)
- `length_bp` — sequence length in base pairs

**Verify:**
- All 5 uploaded files appear in the table
- Each has a pLIN code with 6 numbers separated by dots
- The Inc type is detected (not "Unknown")

### Step E6. Check the Cladogram Tab

Click the **"Cladogram"** tab.

**You should see:**
- A tree diagram (dendrogram) showing relationships between the 5 plasmids
- Branches colored by Inc group
- Interactive features: hover over branches to see distances

**Try:**
- Zoom in/out with your mouse scroll wheel
- Hover over leaf nodes to see plasmid names
- If available, switch between "Rectangular" and "Circular" views

### Step E7. Check the AMR Analysis Tab (if AMRFinderPlus is installed)

Click the **"AMR Analysis"** tab.

**If AMRFinderPlus is installed, you should see:**
- A table of detected AMR genes per plasmid
- Drug class breakdown (pie chart or bar chart)
- Resistance gene heatmap

**If AMRFinderPlus is NOT installed:**
- You will see a message saying the tool is not available
- This is expected — AMR analysis is optional

### Step E8. Check the Epidemiology Tab

Click the **"Epidemiology"** tab.

**You should see:**
- Mobility classification results (if MOBsuite is installed)
- Cluster analysis showing which plasmids group together
- Outbreak detection alerts (if any plasmids share identical pLIN codes)

### Step E9. Export Results

Click the **"Export"** tab.

1. Click **"Download pLIN Assignments (TSV)"**
   - A `.tsv` file downloads to your computer
   - Open it in Excel, Google Sheets, or a text editor
   - **Verify:** it has columns `plasmid_id`, `inc_type`, `length_bp`, `pLIN`, `bin_A` through `bin_F`
   - **Verify:** all 5 plasmids are listed with numeric pLIN codes

2. Click **"Download All Results (ZIP)"** (if available)
   - A `.zip` file downloads
   - Extract it and verify it contains `.tsv` files

### Step E10. Upload Your Own Data (Optional)

If you have your own plasmid FASTA files:

1. Click **"Browse files"** again
2. Upload your FASTA files (`.fasta`, `.fa`, or `.fna` format)
3. Click **"Run Analysis"**
4. Check that results appear without errors

**Supported input:**
- Standard FASTA format (header line starting with `>` followed by DNA sequence)
- One or more sequences per file
- File size up to 200 MB
- Any number of files

---

## F. Verify Each Feature Works

Use this checklist to confirm every feature of pLIN works correctly.

### Core Features (No External Tools Required)

| Step | What to Do | Expected Result | Pass? |
|---|---|---|---|
| F1 | Upload 5 IncX FASTA files | Files accepted, count shown | |
| F2 | Click "Run Analysis" | Analysis completes, no errors | |
| F3 | Check Results tab | Table with pLIN codes for all 5 files | |
| F4 | Check pLIN format | 6 numbers separated by dots (e.g., `3.5.12.45.201.3050`) | |
| F5 | Check Inc type detection | All 5 files classified (not "Unknown") | |
| F6 | Check Cladogram tab | Tree diagram renders with 5 leaves | |
| F7 | Zoom/pan the cladogram | Interactive controls work | |
| F8 | Check Overview tab | Threshold table and metrics displayed | |
| F9 | Download TSV from Export tab | File downloads, opens in Excel/text editor | |
| F10 | Download ZIP from Export tab | ZIP file downloads with multiple files | |
| F11 | Change Inc Group to "IncX3" in sidebar | Results update when re-run | |
| F12 | Change Linkage to "Complete" | Cladogram changes when re-run | |

### Optional Features (Require External Tools)

These only work if the corresponding tool is installed. The sidebar shows which tools are available.

| Step | Feature | Tool Required | What to Check |
|---|---|---|---|
| F13 | AMR gene detection | AMRFinderPlus | AMR Analysis tab shows gene table |
| F14 | Gene prediction | Prodigal | Gene annotations in results |
| F15 | Mobility typing | MOBsuite | Relaxase family + MPF type in Epidemiology tab |
| F16 | Fast ANI estimation | Mash | Distance matrix in Epidemiology tab |
| F17 | Precise ANI calculation | FastANI | ANI concordance validation |
| F18 | SNP sub-typing | minimap2 | SNP counts within clusters |
| F19 | CRISPR host prediction | MinCED + BLAST+ | CRISPR Host tab shows spacer results |
| F20 | Chromosomal typing | mlst | MLST results in Epidemiology tab |

### New Analytical Modules (v3.0)

These features run automatically as part of the analysis.

| Step | Feature | What to Check |
|---|---|---|
| F21 | Assembly completeness | Completeness scores in Quality Assessment section |
| F22 | Database coverage | Novelty flags in Classification Results section |
| F23 | Novel Inc discovery | Novel group clustering (if unknown plasmids found) |
| F24 | Recombination detection | Recombination signals in Epidemiology tab |
| F25 | Evolutionary rate estimation | SNP accumulation rates (requires dated metadata) |
| F26 | Cluster stability | Bootstrap support scores in Classification Results |
| F27 | MGE boundary detection | Gene architecture map (requires Prodigal) |

---

## G. Quick Validation Test (5 Minutes)

If you are short on time, do only these steps:

1. **Start pLIN** (using your chosen method from Sections A-D)
2. **Open** `http://localhost:8501` in your browser
3. **Upload** 3 files from `test_plasmids/IncX/`:
   - `IncX3_JN247852.fasta`
   - `SP12_P2.fasta`
   - `SP42_P2.fasta`
4. **Click** "Run Analysis"
5. **Check** the Results tab: you should see 3 rows with pLIN codes
6. **Check** the Cladogram tab: you should see a tree with 3 leaves
7. **Download** the TSV file from the Export tab
8. **Open** the TSV file — verify it has columns and data

If all 8 steps work without errors, the core tool is validated.

---

## H. Troubleshooting

### "Command not found: conda"

**Cause:** Conda is not in your PATH.

**Fix (macOS):**

```bash
source ~/miniconda3/bin/activate
conda init zsh
```

Close and reopen your terminal.

**Fix (Linux):**

```bash
source ~/miniconda3/bin/activate
conda init bash
```

Close and reopen your terminal.

**Fix (Windows):** Open **Anaconda Prompt** from the Start Menu instead of Command Prompt or PowerShell.

---

### "ModuleNotFoundError: No module named 'streamlit'"

**Cause:** Python packages are not installed, or you are not in the correct environment.

**Fix:**

```bash
conda activate pLIN_tools
pip install -r requirements.txt
```

---

### "Port 8501 is already in use"

**Cause:** Another application is using port 8501.

**Fix:** Use a different port:

```bash
streamlit run plin_app.py --server.port 8502
```

Then open `http://localhost:8502` in your browser.

---

### "IndexError" or "FileNotFoundError" When Running Analysis

**Cause:** Data files may be missing or corrupted.

**Fix:** Verify these files exist:

```bash
ls data/inc_classifier.npz
ls data/inc_centroids.npz
ls output/pLIN_assignments.tsv
```

If any file is missing, re-clone the repository:

```bash
git clone https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification.git
```

---

### Browser Shows Blank Page

**Cause:** Streamlit is still starting up.

**Fix:**
1. Wait 10 seconds and refresh the page
2. Check the terminal — it should say "You can now view your Streamlit app in your browser"
3. Try `http://127.0.0.1:8501` instead of `http://localhost:8501`

---

### Docker Build Fails

**Cause:** Docker daemon is not running, or insufficient disk space.

**Fix:**
1. Make sure Docker Desktop is running (look for the whale icon in your system tray)
2. Check disk space: you need at least 10 GB free
3. Try rebuilding with no cache:
   ```bash
   docker build --no-cache --build-arg INSTALL_BIOTOOLS=true -t plin-full .
   ```

---

### "Permission denied" on Linux

**Cause:** Script is not executable.

**Fix:**

```bash
chmod +x pLIN_Linux.sh
./pLIN_Linux.sh
```

---

### AMRFinderPlus / Prodigal / Other Tool "Not Found"

**Cause:** Optional bioinformatics tools are not installed. This is normal.

**What happens:** The app works without these tools. Features that need them will show a message saying the tool is not available. Core pLIN classification always works.

**To install missing tools:**

```bash
conda activate pLIN_tools
conda install -c bioconda -c conda-forge ncbi-amrfinderplus mash fastani minimap2 minced blast prodigal mlst -y
pip install mob_suite
amrfinder --update
```

---

### Windows: Bioinformatics Tools Not Working

**Cause:** Most bioinformatics tools are designed for Unix (macOS/Linux).

**Fix:** Use one of these options:
1. **WSL2** (recommended): Run `wsl --install` in PowerShell, then use Linux instructions inside WSL
2. **Docker** (easiest): Follow Section D
3. **Native Windows:** Core pLIN classification works, but AMR detection and other advanced features require WSL2 or Docker

---

## I. Environment Check Command

To see a summary of your setup, run:

**macOS:**

```bash
./pLIN_macOS.sh --check
```

**Linux:**

```bash
./pLIN_Linux.sh --check
```

**Windows:**

```cmd
pLIN_Windows.bat --check
```

This shows:
- Your Python version
- Which data files are present
- Which Python packages are installed
- Which bioinformatics tools are available

---

## J. Stopping and Restarting

### Stop pLIN

- **Native install:** Press **Ctrl+C** in the terminal where pLIN is running
- **Docker Compose:** Run `docker compose down`
- **Docker manual:** Run `docker stop plin-app`

### Restart pLIN Later

- **macOS:** Open Terminal, navigate to the pLIN folder, run `./pLIN_macOS.sh --launch`
- **Linux:** Open Terminal, navigate to the pLIN folder, run `./pLIN_Linux.sh --launch`
- **Windows:** Open Anaconda Prompt, navigate to the pLIN folder, run `pLIN_Windows.bat --launch`
- **Docker Compose:** Run `docker compose up`

### Uninstall pLIN

- **macOS:** `./pLIN_macOS.sh --uninstall`
- **Linux:** `./pLIN_Linux.sh --uninstall`
- **Windows:** `pLIN_Windows.bat --uninstall`
- **Docker:** `docker rmi plin-full` and `docker rmi plin:2.1.0`

This removes the conda environment and virtual environment. Your data files are preserved.

---

## K. File Format Reference

### Input: FASTA Format

pLIN accepts standard FASTA files. Each file should look like this:

```
>plasmid_name some optional description
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
```

- Line 1: starts with `>` followed by the sequence name
- Remaining lines: DNA sequence (A, T, C, G characters)
- File extension: `.fasta`, `.fa`, or `.fna`
- One or more sequences per file
- Maximum upload size: 200 MB

### Input: Metadata (Optional)

A CSV or TSV file with epidemiological data:

```
plasmid_id,collection_date,location,patient_id,source
SP12_P2,2023-03-15,Hospital_A,Patient_001,blood
SP42_P2,2023-03-18,Hospital_A,Patient_002,urine
SP30_P3,2023-04-01,Hospital_B,Patient_003,wound
```

- Required column: `plasmid_id` (must match uploaded FASTA file names)
- Optional columns: `collection_date`, `location`, `patient_id`, `source`
- Dates should be in `YYYY-MM-DD` format

### Output: pLIN Assignments TSV

The main output file has these columns:

| Column | Description | Example |
|---|---|---|
| `plasmid_id` | Plasmid name from FASTA header | `SP12_P2` |
| `inc_type` | Detected incompatibility group | `IncX3` |
| `length_bp` | Sequence length in base pairs | `46894` |
| `pLIN` | Full 6-level classification code | `3.5.12.45.201.3050` |
| `bin_A` | Level 1 cluster (broadest) | `3` |
| `bin_B` | Level 2 cluster | `5` |
| `bin_C` | Level 3 cluster | `12` |
| `bin_D` | Level 4 cluster | `45` |
| `bin_E` | Level 5 cluster | `201` |
| `bin_F` | Level 6 cluster (finest / strain-level) | `3050` |

---

## L. System Requirements Summary

| Component | Minimum | Recommended |
|---|---|---|
| **Operating System** | macOS 13+, Ubuntu 20.04+, Windows 10 | macOS 14+, Ubuntu 22.04+, Windows 11 |
| **Python** | 3.10 | 3.11 or 3.12 |
| **RAM** | 4 GB | 8 GB |
| **Disk Space** | 5 GB | 10 GB |
| **Browser** | Chrome, Firefox, Safari, Edge | Chrome or Firefox (latest) |
| **Internet** | Required for first-time setup | Not needed after setup |

---

## M. Contact and Citation

**Repository:** https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification

**License:** GPL-3.0 with mandatory citation clause

**Required Citation:**

> Xavier, B. (2025). pLIN: A Plasmid Lineage Identification Number System for Hierarchical, Permanent Classification of Bacterial Plasmids Integrated with Antimicrobial Resistance Gene Surveillance. https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification
