# pLIN Quick Start Guide — Linux

**For reviewers with no bioinformatics experience.**
**Time needed: 20-30 minutes (mostly waiting for downloads).**

---

## What is pLIN?

pLIN is a web application that runs on your computer. You will open it in your web browser (Chrome, Firefox, or Edge) — just like any website, except it runs locally on your machine. Nothing is uploaded to the internet.

You give it plasmid DNA sequence files, and it:
- Classifies them into groups (Inc/Rep types)
- Assigns a 6-level code (like `3.5.12.45.201.3050`)
- Detects antimicrobial resistance (AMR) genes
- Shows interactive trees and charts

---

## Before You Start

You need:
- **Ubuntu 20.04+**, **Debian 11+**, **Fedora 36+**, **CentOS 8+**, **Arch Linux**, or **RHEL 8+**
- At least **5 GB** of free disk space
- **sudo access** (administrator privileges — you will be asked for your password)
- An **internet connection** (only needed during setup)
- About **20-30 minutes** of time

---

## Step 1: Open Terminal

- **Ubuntu/Debian:** Press **Ctrl + Alt + T**
- **Fedora/RHEL:** Click "Activities" (top-left), type **Terminal**, click it
- **Any Linux:** Right-click the desktop and choose **"Open Terminal"**

A window with a blinking cursor will appear. This is your Terminal.

> **Tip:** You will copy-paste commands from this guide. To paste in Terminal, press **Ctrl + Shift + V** (not Ctrl + V).

---

## Step 2: Install Basic Tools

First, install git, curl, and wget (you may already have them).

**Ubuntu / Debian:**
```bash
sudo apt update && sudo apt install -y git curl wget unzip
```

**Fedora / RHEL:**
```bash
sudo dnf install -y git curl wget unzip
```

**Arch Linux:**
```bash
sudo pacman -Sy --noconfirm git curl wget unzip
```

When asked for your password, type it and press Enter. **Nothing appears on screen while you type your password** — this is normal.

---

## Step 3: Unzip the pLIN Package

You received a file called `pLIN_v3.0.0_Linux.zip`.

1. Find the ZIP file (likely in your **Downloads** folder)
2. Unzip it using Terminal:

```bash
cd ~/Downloads
unzip pLIN_v3.0.0_Linux.zip -d pLIN_v3.0.0_Linux
cd pLIN_v3.0.0_Linux
```

> **Note:** If you saved the ZIP somewhere else, adjust the path accordingly.

3. Verify you are in the right folder:

```bash
ls plin_app.py
```

You should see `plin_app.py` printed. If you see "No such file or directory", you are in the wrong folder.

---

## Step 4: Install Miniconda (Python Environment Manager)

Miniconda manages Python and scientific software packages.

**Check if you already have it:**

```bash
conda --version
```

- If you see `conda 24.x.x` — **skip to Step 5**
- If you see "command not found" — install it:

**For 64-bit Intel/AMD systems (most common):**
```bash
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -p $HOME/miniconda3
rm /tmp/miniconda.sh
```

**For ARM64 systems (Raspberry Pi, some servers):**
```bash
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -o /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -p $HOME/miniconda3
rm /tmp/miniconda.sh
```

> **Not sure which one?** Run `uname -m`. If it says `x86_64`, use the first option. If it says `aarch64`, use the second.

Now initialize conda:

```bash
~/miniconda3/bin/conda init bash
```

**Important: Close your Terminal window completely and open a new one.**

Navigate back to the pLIN folder:

```bash
cd ~/Downloads/pLIN_v3.0.0_Linux
```

**Verify it worked:**

```bash
conda --version
```

You should see `conda 24.x.x` or similar.

---

## Step 5: Run the Installer

This is the main installation step. It will download and install everything pLIN needs.

```bash
chmod +x pLIN_Linux.sh
./pLIN_Linux.sh --install
```

**What to expect:**
- It may ask for your **sudo password** (to install system packages like Java)
- You will see `[INFO]` and `[OK]` messages as things install
- Some items may show `[WARN] Could not install` — **this is fine**, these are optional tools
- **This takes 15-30 minutes** — let it run, do not close Terminal
- When finished, you will see:

```
  [OK]   Installation complete!
```

> **If you see errors about "command not found: conda":**
> Run `source ~/miniconda3/bin/activate` and try the install command again.

---

## Step 6: Launch pLIN

```bash
./pLIN_Linux.sh --launch
```

**What happens:**
1. Terminal shows "pLIN is starting at: http://localhost:8501"
2. Your web browser opens automatically
3. You see the **pLIN Classifier** web interface

> **If the browser does not open automatically:**
> Open Firefox or Chrome and go to: **http://localhost:8501**

> **If you see a blank page:** Wait 10 seconds and refresh the page (Ctrl + R or F5).

---

## Step 7: Test with Sample Data

The package includes 22 test plasmid files. Let's use a few to verify everything works.

### 7a. Upload Test Files

1. In the pLIN web interface, click **"Browse files"** (in the upload area)
2. Navigate to the `test_plasmids/IncX/` folder inside your pLIN folder
   - Typically: `Downloads > pLIN_v3.0.0_Linux > test_plasmids > IncX`
3. Select these 3 files:
   - `IncX3_JN247852.fasta`
   - `SP12_P2.fasta`
   - `SP42_P2.fasta`
4. Click **Open**

You should see "3 files uploaded" with the file names listed.

### 7b. Run the Analysis

1. Leave all sidebar settings at their defaults
2. Click the **"Run Analysis"** button
3. Wait 10-30 seconds — you will see a progress bar

### 7c. Check the Results

**Results tab:**
- You should see a table with 3 rows
- Each row has a `pLIN` code (6 numbers separated by dots, like `3.5.12.45.201.3050`)
- The `inc_type` column should show an IncX variant (not "Unknown")

**Cladogram tab:**
- You should see a tree diagram with 3 branches
- You can zoom in/out with your mouse scroll wheel
- Hover over branches to see distances

**Export tab:**
- Click **"Download pLIN Assignments (TSV)"**
- A `.tsv` file downloads — open it in a text editor or LibreOffice Calc
- Verify it has columns: `plasmid_id`, `inc_type`, `length_bp`, `pLIN`

---

## Step 8: Stop pLIN

When you are done testing:

1. Go back to the Terminal window where pLIN is running
2. Press **Ctrl + C** (hold Control and press C)
3. Terminal returns to the normal command prompt

---

## Restarting pLIN Later

Next time you want to use pLIN, you only need two commands:

```bash
cd ~/Downloads/pLIN_v3.0.0_Linux
./pLIN_Linux.sh --launch
```

No reinstallation needed.

---

## Troubleshooting

### "command not found: conda"

Run this, then try again:
```bash
source ~/miniconda3/bin/activate
conda init bash
```
Close and reopen Terminal.

### "Permission denied" when running the script

```bash
chmod +x pLIN_Linux.sh
./pLIN_Linux.sh --install
```

### "ModuleNotFoundError: No module named 'streamlit'"

```bash
conda activate pLIN_tools
pip install -r requirements.txt
```

### "Port 8501 is already in use"

Another application is using that port. Use a different one:
```bash
conda activate pLIN_tools
streamlit run plin_app.py --server.port 8502
```
Then open **http://localhost:8502** in your browser.

### Browser shows a blank page

1. Wait 10 seconds, then refresh (Ctrl + R or F5)
2. Try **http://127.0.0.1:8501** instead
3. Check Terminal for error messages

### "E: Unable to locate package" (Ubuntu/Debian)

Your package list is outdated. Run:
```bash
sudo apt update
```
Then try the command again.

### AMR Analysis tab says "tool not available"

This means AMRFinderPlus did not install (it is optional). Core pLIN classification works perfectly without it. To install manually:
```bash
conda activate pLIN_tools
conda install -c bioconda -c conda-forge ncbi-amrfinderplus -y
amrfinder --update
```

---

## What Each File Does (Reference)

| File | Purpose |
|------|---------|
| `plin_app.py` | The main application (do not edit) |
| `pLIN_Linux.sh` | Installer and launcher script for Linux |
| `setup_pLIN.py` | Cross-platform installer (alternative) |
| `requirements.txt` | List of Python packages needed |
| `data/inc_classifier.npz` | Classification model (trained on 8,077 plasmids) |
| `data/inc_centroids.npz` | Group centroids for distance calculation |
| `output/pLIN_assignments.tsv` | Reference database (79,305 plasmid assignments) |
| `test_plasmids/IncX/*.fasta` | 22 sample plasmid sequences for testing |
| `REVIEWER_GUIDE.md` | Detailed feature validation checklist |
| `Dockerfile` | For Docker-based installation (alternative) |

---

## Alternative: Docker Installation (Easiest)

If you prefer Docker, this avoids all manual dependency installation:

1. Install Docker:
   - **Ubuntu/Debian:** `sudo apt update && sudo apt install -y docker.io docker-compose-plugin`
   - **Fedora/RHEL:** `sudo dnf install -y docker docker-compose-plugin && sudo systemctl start docker`
   - Then: `sudo usermod -aG docker $USER` and **log out and log back in**

2. Run pLIN:
   ```bash
   cd ~/Downloads/pLIN_v3.0.0_Linux
   docker compose up --build
   ```

3. Open **http://localhost:8501** in your browser

4. Stop with: `docker compose down`

First build takes 10-20 minutes.
