 # pLIN Quick Start Guide — Windows
    
**For reviewers with no bioinformatics experience.**
**Time needed: 25-40 minutes (mostly waiting for downloads).**

---

## What is pLIN?

pLIN is a web application that runs on your computer. You will open it in your web browser (Chrome, Edge, or Firefox) — just like any website, except it runs locally on your PC. Nothing is uploaded to the internet.

You give it plasmid DNA sequence files, and it:
- Classifies them into groups (Inc/Rep types)
- Assigns a 6-level code (like `3.5.12.45.201.3050`)
- Detects antimicrobial resistance (AMR) genes
- Shows interactive trees and charts

---

## Before You Start

You need:
- **Windows 10** (version 1903 or newer) or **Windows 11**
- At least **5 GB** of free disk space
- An **internet connection** (only needed during setup)
- About **25-40 minutes** of time

---

## Step 1: Install Git for Windows

Git is a tool for downloading software. You may already have it.

**Check if you already have it:**
1. Press **Windows key + R**, type `cmd`, press Enter
2. Type `git --version` and press Enter
3. If you see `git version 2.x.x` — **skip to Step 2**

**If you don't have it:**
1. Open your web browser and go to: **https://git-scm.com/download/win**
2. The download should start automatically
3. Run the downloaded installer
4. Click **Next** on every screen (the default settings are fine)
5. Click **Install**
6. Click **Finish**

---

## Step 2: Install Miniconda

Miniconda manages Python and scientific software. You may already have it.

**Check if you already have it:**
1. Click the **Start Menu** (Windows icon at the bottom-left)
2. Type **Anaconda Prompt**
3. If you see "Anaconda Prompt" in the results — **skip to Step 3**

**If you don't have it:**
1. Open your web browser and go to: **https://docs.conda.io/en/latest/miniconda.html**
2. Download **Miniconda3 Windows 64-bit** (the `.exe` file)
3. Run the downloaded installer
4. Click **Next**, then **I Agree**
5. Choose **Just Me** and click **Next**
6. Keep the default install location and click **Next**
7. On the **Advanced Options** screen:
   - **CHECK the box** that says **"Add Miniconda3 to my PATH environment variable"**
   - (It says "Not recommended" — ignore that warning, we need this)
8. Click **Install**
9. Wait for installation to finish
10. Click **Finish**

---

## Step 3: Open Anaconda Prompt

This is the command window you will use for all remaining steps.

1. Click the **Start Menu** (Windows icon)
2. Type **Anaconda Prompt**
3. Click **Anaconda Prompt** (the one with the green icon)

A black window with a blinking cursor will appear. This is your command prompt.

> **Important:** Always use **Anaconda Prompt** for pLIN commands — not regular Command Prompt or PowerShell.

> **Tip:** You will copy-paste commands from this guide. To paste in Anaconda Prompt, **right-click** inside the window.

---

## Step 4: Unzip the pLIN Package

You received a file called `pLIN_v3.0.0_Windows.zip`.

1. Find the ZIP file (likely in your **Downloads** folder)
2. **Right-click** the ZIP file and choose **"Extract All..."**
3. Click **Extract** (this creates a folder called `pLIN_v3.0.0_Windows`)
4. In Anaconda Prompt, navigate to that folder:

```cmd
cd %USERPROFILE%\Downloads\pLIN_v3.0.0_Windows
```

> **Note:** If you saved the ZIP somewhere else, adjust the path. For example, if it is on your Desktop:
> ```cmd
> cd %USERPROFILE%\Desktop\pLIN_v3.0.0_Windows
> ```

5. Verify you are in the right folder:

```cmd
dir plin_app.py
```

You should see `plin_app.py` listed. If you see "File Not Found", you are in the wrong folder.

---

## Step 5: Run the Installer

This is the main installation step. It will download and install everything pLIN needs.

```cmd
pLIN_Windows.bat --install
```

**What to expect:**
- You will see `[INFO]` and `[OK]` messages as things install
- Some items may show `[WARN] Could not install` — **this is fine**, these are optional tools
- **This takes 15-30 minutes** — let it run, do not close the window
- When finished, you will see:

```
  Installation Complete!
```

Then press any key when prompted.

> **Note about Windows:** Some bioinformatics tools (AMRFinderPlus, Prodigal) have limited Windows support. The core pLIN classification works perfectly, but for full AMR analysis, see the "Full Feature Support" section at the bottom.

---

## Step 6: Launch pLIN

```cmd
pLIN_Windows.bat --launch
```

**What happens:**
1. The prompt shows "pLIN is starting at: http://localhost:8501"
2. Your web browser opens automatically
3. You see the **pLIN Classifier** web interface

> **If the browser does not open automatically:**
> Open Chrome, Edge, or Firefox and go to: **http://localhost:8501**

> **If you see a blank page:** Wait 10 seconds and refresh the page (F5 key or Ctrl + R).

---

## Step 7: Test with Sample Data

The package includes 22 test plasmid files. Let's use a few to verify everything works.

### 7a. Upload Test Files

1. In the pLIN web interface, click **"Browse files"** (in the upload area)
2. Navigate to the `test_plasmids\IncX\` folder inside your pLIN folder
   - Typically: `Downloads > pLIN_v3.0.0_Windows > test_plasmids > IncX`
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
- A `.tsv` file downloads — open it in Excel or Notepad
- Verify it has columns: `plasmid_id`, `inc_type`, `length_bp`, `pLIN`

---

## Step 8: Stop pLIN

When you are done testing:

1. Go back to the Anaconda Prompt window where pLIN is running
2. Press **Ctrl + C** (hold Control and press C)
3. If it asks "Terminate batch job? (Y/N)", type **Y** and press Enter
4. The prompt returns to normal

---

## Restarting pLIN Later

Next time you want to use pLIN:

1. Open **Anaconda Prompt** from the Start Menu
2. Navigate to the pLIN folder:

```cmd
cd %USERPROFILE%\Downloads\pLIN_v3.0.0_Windows
```

3. Launch:

```cmd
pLIN_Windows.bat --launch
```

No reinstallation needed.

---

## Troubleshooting

### "conda is not recognized"

Miniconda was not added to PATH. Two options:

**Option A:** Use Anaconda Prompt (not regular Command Prompt)

**Option B:** Reinstall Miniconda and CHECK the "Add to PATH" box this time.

### "'pLIN_Windows.bat' is not recognized"

You are not in the correct folder. Run:
```cmd
cd %USERPROFILE%\Downloads\pLIN_v3.0.0_Windows
dir plin_app.py
```
If you see `plin_app.py`, try the command again.

### "ModuleNotFoundError: No module named 'streamlit'"

```cmd
conda activate pLIN_tools
pip install -r requirements.txt
```

### "Port 8501 is already in use"

Another application is using that port. Use a different one:
```cmd
conda activate pLIN_tools
streamlit run plin_app.py --server.port 8502
```
Then open **http://localhost:8502** in your browser.

### Browser shows a blank page

1. Wait 10 seconds, then refresh (F5)
2. Try **http://127.0.0.1:8501** instead
3. Check Anaconda Prompt for error messages

### AMR Analysis tab says "tool not available"

This is expected on native Windows. Most bioinformatics tools are designed for macOS/Linux. Core pLIN classification works perfectly without them.

For full feature support, see the next section.

---

## Full Feature Support on Windows (Optional)

For all bioinformatics tools (AMR detection, CRISPR analysis, etc.), you have two options:

### Option A: Windows Subsystem for Linux (WSL2)

WSL2 lets you run Linux inside Windows. This gives you full tool support.

1. Open **PowerShell as Administrator**:
   - Right-click the Start Menu > **Terminal (Admin)** or **Windows PowerShell (Admin)**
2. Run:
   ```powershell
   wsl --install
   ```
3. **Restart your computer**
4. After restart, Ubuntu will open automatically and ask you to create a username and password
5. In the Ubuntu window, run:
   ```bash
   cd /mnt/c/Users/YOUR_USERNAME/Downloads/pLIN_v3.0.0_Windows
   chmod +x pLIN_Linux.sh
   ./pLIN_Linux.sh --install
   ./pLIN_Linux.sh --launch
   ```
   (Replace `YOUR_USERNAME` with your actual Windows username)

### Option B: Docker Desktop

1. Download Docker Desktop from: **https://www.docker.com/products/docker-desktop/**
2. Run the installer — when asked, enable **WSL2 backend**
3. **Restart your computer**
4. Open Docker Desktop (wait until it says "Docker is running")
5. Open Anaconda Prompt and run:
   ```cmd
   cd %USERPROFILE%\Downloads\pLIN_v3.0.0_Windows
   docker compose up --build
   ```
6. Open **http://localhost:8501** in your browser
7. Stop with: `docker compose down`

---

## What Each File Does (Reference)

| File | Purpose |
|------|---------|
| `plin_app.py` | The main application (do not edit) |
| `pLIN_Windows.bat` | Installer and launcher script for Windows |
| `setup_pLIN.py` | Cross-platform installer (alternative) |
| `requirements.txt` | List of Python packages needed |
| `data\inc_classifier.npz` | Classification model (trained on 8,077 plasmids) |
| `data\inc_centroids.npz` | Group centroids for distance calculation |
| `output\pLIN_assignments.tsv` | Reference database (79,305 plasmid assignments) |
| `test_plasmids\IncX\*.fasta` | 22 sample plasmid sequences for testing |
| `REVIEWER_GUIDE.md` | Detailed feature validation checklist |
| `Dockerfile` | For Docker-based installation (alternative) |
