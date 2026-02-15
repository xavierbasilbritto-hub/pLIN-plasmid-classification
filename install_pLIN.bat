@echo off
REM ═══════════════════════════════════════════════════════════════════════════
REM pLIN Tool — Installation Script (Windows)
REM Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
REM ═══════════════════════════════════════════════════════════════════════════

echo ═══════════════════════════════════════════════════════════════
echo  pLIN Tool — Installer for Windows
echo ═══════════════════════════════════════════════════════════════
echo.

set ENV_NAME=pLIN_tools

REM ── Check for conda ──────────────────────────────────────────────────────
where conda >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo [!] Conda not found.
    echo     Please install Miniconda from:
    echo     https://docs.conda.io/en/latest/miniconda.html
    echo.
    echo     During installation, check "Add Miniconda3 to my PATH"
    echo     Then restart this script from Anaconda Prompt.
    echo.
    pause
    exit /b 1
)

echo [OK] Conda found.

REM ── Create / activate conda environment ──────────────────────────────────
echo.
echo [1/5] Setting up conda environment: %ENV_NAME%

conda env list | findstr /C:"%ENV_NAME%" >nul 2>&1
if %ERRORLEVEL% equ 0 (
    echo       Environment '%ENV_NAME%' already exists. Activating...
) else (
    echo       Creating new environment with Python 3.11...
    conda create -n %ENV_NAME% python=3.11 -y
)

call conda activate %ENV_NAME%
echo       Python:
python --version

REM ── Install Python dependencies ──────────────────────────────────────────
echo.
echo [2/5] Installing Python dependencies...
pip install --upgrade pip
pip install streamlit numpy pandas scipy biopython scikit-learn matplotlib seaborn plotly python-pptx requests
echo       Core Python packages installed.

REM ── Install bioinformatics tools via conda ───────────────────────────────
echo.
echo [3/5] Installing bioinformatics tools...

conda config --add channels bioconda 2>nul
conda config --add channels conda-forge 2>nul

echo       Installing AMRFinderPlus...
conda install -c bioconda -c conda-forge ncbi-amrfinderplus -y 2>nul
if %ERRORLEVEL% equ 0 (echo       [OK] AMRFinderPlus installed) else (echo       [SKIP] AMRFinderPlus - optional)

echo       Installing Mash...
conda install -c bioconda mash -y 2>nul
if %ERRORLEVEL% equ 0 (echo       [OK] Mash installed) else (echo       [SKIP] Mash - optional)

echo       Installing FastANI...
conda install -c bioconda fastani -y 2>nul
if %ERRORLEVEL% equ 0 (echo       [OK] FastANI installed) else (echo       [SKIP] FastANI - optional)

echo       Installing minimap2...
conda install -c bioconda minimap2 -y 2>nul
if %ERRORLEVEL% equ 0 (echo       [OK] minimap2 installed) else (echo       [SKIP] minimap2 - optional)

echo       Installing MinCED...
conda install -c bioconda minced -y 2>nul
if %ERRORLEVEL% equ 0 (echo       [OK] MinCED installed) else (echo       [SKIP] MinCED - optional)

echo       Installing BLAST+...
conda install -c bioconda blast -y 2>nul
if %ERRORLEVEL% equ 0 (echo       [OK] BLAST+ installed) else (echo       [SKIP] BLAST+ - optional)

echo       Installing Prodigal...
conda install -c bioconda prodigal -y 2>nul
if %ERRORLEVEL% equ 0 (echo       [OK] Prodigal installed) else (echo       [SKIP] Prodigal - optional)

REM ── Install MOBsuite via pip ─────────────────────────────────────────────
echo.
echo [4/5] Installing MOBsuite...
pip install mob_suite 2>nul
if %ERRORLEVEL% equ 0 (echo       [OK] MOBsuite installed) else (echo       [SKIP] MOBsuite - optional)

REM ── Update AMRFinderPlus database ────────────────────────────────────────
echo.
echo [5/5] Updating AMRFinderPlus database...
amrfinder --update 2>nul
if %ERRORLEVEL% equ 0 (echo       [OK] Database updated) else (echo       [SKIP] Run manually: amrfinder --update)

REM ── Summary ──────────────────────────────────────────────────────────────
echo.
echo ═══════════════════════════════════════════════════════════════
echo  Installation Complete!
echo ═══════════════════════════════════════════════════════════════
echo.
echo  To launch pLIN:
echo    conda activate %ENV_NAME%
echo    streamlit run plin_app.py
echo.
echo  The app will open at: http://localhost:8501
echo.
echo  Optional: Install Ollama for AI chatbot:
echo    Download from: https://ollama.com/download/windows
echo    Then run: ollama pull llama3.2
echo.
echo ═══════════════════════════════════════════════════════════════
pause
