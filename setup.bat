@echo off
REM ============================================================================
REM pLIN Setup Script — Windows
REM Creates a Python virtual environment and installs all dependencies.
REM ============================================================================

echo ============================================================
echo   pLIN — Plasmid Life Identification Number System
echo   Setup Script (Windows)
echo ============================================================
echo.

REM ── 1. Python virtual environment ──────────────────────────────────────────
echo [1/3] Setting up Python virtual environment ...

if exist .venv (
    echo   .venv already exists. Skipping creation.
) else (
    python -m venv .venv
    echo   Created .venv
)

call .venv\Scripts\activate.bat
echo   Activated .venv

echo.
echo [2/3] Installing Python dependencies ...
pip install --upgrade pip -q
pip install -r requirements.txt -q
echo   All Python packages installed.

REM ── 2. Create directory structure ──────────────────────────────────────────
echo.
echo [3/3] Creating directory structure ...
if not exist "plasmid_sequences_for_training\IncFII\fastas" mkdir "plasmid_sequences_for_training\IncFII\fastas"
if not exist "plasmid_sequences_for_training\IncN\fastas" mkdir "plasmid_sequences_for_training\IncN\fastas"
if not exist "plasmid_sequences_for_training\IncX1\fastas" mkdir "plasmid_sequences_for_training\IncX1\fastas"
if not exist "output\amrfinder" mkdir "output\amrfinder"
if not exist "output\integrated" mkdir "output\integrated"
if not exist "output\figures" mkdir "output\figures"
echo   Directory structure ready.

echo.
echo ============================================================
echo   Setup Complete!
echo ============================================================
echo.
echo   Next steps:
echo   1. Place your FASTA files in:
echo        plasmid_sequences_for_training\IncFII\fastas\
echo        plasmid_sequences_for_training\IncN\fastas\
echo        plasmid_sequences_for_training\IncX1\fastas\
echo.
echo   2. Run the full pipeline:
echo        run_all.bat
echo.
echo   NOTE: AMRFinderPlus is not natively available on Windows.
echo   Use WSL (Windows Subsystem for Linux) for AMRFinderPlus:
echo     wsl conda install -c bioconda ncbi-amrfinderplus
echo.
echo ============================================================
pause
