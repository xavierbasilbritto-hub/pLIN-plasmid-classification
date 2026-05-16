# pLIN Quick Start Guide — macOS

**For reviewers with no bioinformatics experience.**
**Time needed: 20-30 minutes (mostly waiting for downloads).**

---

## What is pLIN?

pLIN is a web application that runs on your computer. You will open it in your web browser (Safari, Chrome, or Firefox) — just like any website, except it runs locally on your Mac. Nothing is uploaded to the internet.

You give it plasmid DNA sequence files, and it:
- Classifies them into groups (Inc/Rep types)
- Assigns a 6-level code (like `3.5.12.45.201.3050`)
- Detects antimicrobial resistance (AMR) genes
- Shows interactive trees and charts

---

## Before You Start

You need:
- A Mac running **macOS 13 (Ventura) or newer** (check: Apple menu > About This Mac)
- At least **5 GB** of free disk space
- An **internet connection** (only needed during setup)
- About **20-30 minutes** of time

---

## Step 1: Open Terminal

Terminal is a program already installed on every Mac. It lets you type commands.

1. Press **Cmd + Space** on your keyboard (this opens Spotlight search)
2. Type **Terminal**
3. Press **Enter**

A window with a black or white background will appear with a blinking cursor. This is your Terminal.

> **Tip:** You will copy-paste commands from this guide into Terminal. To paste in Terminal, press **Cmd + V**.

---

## Step 2: Unzip the pLIN Package

You received a file called `pLIN_v3.0.0_macOS.zip`.

1. Find the ZIP file (likely in your **Downloads** folder)
2. **Double-click** it to unzip — this creates a folder called `pLIN_v3.0.0_macOS`
3. In Terminal, type this command to navigate into that folder:

```bash
cd ~/Downloads/pLIN_v3.0.0_macOS
```

> **Note:** If you saved the ZIP somewhere else, adjust the path. For example, if it is on your Desktop:
> ```bash
> cd ~/Desktop/pLIN_v3.0.0_macOS
> ```

4. Verify you are in the right folder:

```bash
ls plin_app.py
```

You should see `plin_app.py` printed. If you see "No such file or directory", you are in the wrong folder.

---

## Step 3: Install Homebrew (Package Manager)

Homebrew is a free tool that helps install software on Mac. You may already have it.

**Check if you already have it:**

```bash
brew --version
```

- If you see `Homebrew 4.x.x` — **skip to Step 4**
- If you see "command not found" — install it:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

**What to expect:**
- It will ask for your **Mac password** — type it and press Enter
- **Nothing appears on screen while you type your password** — this is normal, just keep typing
- Wait 2-5 minutes until you see "Installation successful"
- If it tells you to run extra commands (about adding Homebrew to your PATH), **run those commands too**

**Verify it worked:**

```bash
brew --version
```

You should now see `Homebrew 4.x.x`.

---

## Step 4: Install Miniconda (Python Environment Manager)

Miniconda manages Python and scientific software packages.

**Check if you already have it:**

```bash
conda --version
```

- If you see `conda 24.x.x` — **skip to Step 5**
- If you see "command not found" — install it:

```bash
brew install --cask miniconda
```

Wait for it to finish, then run:

```bash
conda init zsh
```

**Important: Close your Terminal window completely and open a new one** (Cmd + Q to quit Terminal, then reopen it from Spotlight).

Navigate back to the pLIN folder:

```bash
cd ~/Downloads/pLIN_v3.0.0_macOS
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
chmod +x pLIN_macOS.sh
./pLIN_macOS.sh --install
```

**What to expect:**
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
./pLIN_macOS.sh --launch
```

**What happens:**
1. Terminal shows "pLIN is starting at: http://localhost:8501"
2. Your web browser opens automatically
3. You see the **pLIN Classifier** web interface

> **If the browser does not open automatically:**
> Open Safari, Chrome, or Firefox and go to: **http://localhost:8501**

> **If you see a blank page:** Wait 10 seconds and refresh the page (Cmd + R).

---

## Step 7: Test with Sample Data

The package includes 22 test plasmid files. Let's use a few to verify everything works.

### 7a. Upload Test Files

1. In the pLIN web interface, click **"Browse files"** (in the upload area)
2. Navigate to the `test_plasmids/IncX/` folder inside your pLIN folder
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
- A `.tsv` file downloads — open it in Excel or any text editor
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
cd ~/Downloads/pLIN_v3.0.0_macOS
./pLIN_macOS.sh --launch
```

No reinstallation needed.

---

## Troubleshooting

### "command not found: conda"

Run this, then try again:
```bash
source ~/miniconda3/bin/activate
conda init zsh
```
Close and reopen Terminal.

### "ModuleNotFoundError: No module named 'streamlit'"

```bash
conda activate pLIN_tools
pip install -r requirements.txt
```

### "Port 8501 is already in use"

Another application is using that port. Use a different one:
```bash
streamlit run plin_app.py --server.port 8502
```
Then open **http://localhost:8502** in your browser.

### Browser shows a blank page

1. Wait 10 seconds, then refresh (Cmd + R)
2. Try **http://127.0.0.1:8501** instead
3. Check Terminal for error messages

### AMR Analysis tab says "tool not available"

This means AMRFinderPlus did not install (it is optional). Core pLIN classification still works perfectly. To install it manually:
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
| `pLIN_macOS.sh` | Installer and launcher script for macOS |
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

If you have Docker Desktop installed, you can skip all the above and just run:

```bash
cd ~/Downloads/pLIN_v3.0.0_macOS
docker compose up --build
```

Then open **http://localhost:8501** in your browser. First build takes 10-20 minutes.

Stop with: `docker compose down`
