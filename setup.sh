#!/bin/bash
# ============================================================================
# pLIN Setup Script — macOS / Linux
# Creates a Python virtual environment and installs all dependencies.
# Optionally installs AMRFinderPlus via conda/mamba.
# ============================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "============================================================"
echo "  pLIN — Plasmid Life Identification Number System"
echo "  Setup Script (macOS / Linux)"
echo "============================================================"
echo ""

# ── 1. Python virtual environment ──────────────────────────────────────────
echo "[1/3] Setting up Python virtual environment ..."

if [ -d ".venv" ]; then
    echo "  .venv already exists. Skipping creation."
else
    python3 -m venv .venv
    echo "  Created .venv"
fi

source .venv/bin/activate
echo "  Activated .venv (Python: $(python3 --version))"

echo ""
echo "[2/3] Installing Python dependencies ..."
pip install --upgrade pip -q
pip install -r requirements.txt -q
echo "  All Python packages installed."

# ── 2. Create directory structure ──────────────────────────────────────────
echo ""
echo "[3/3] Creating directory structure ..."
mkdir -p plasmid_sequences_for_training/IncFII/fastas
mkdir -p plasmid_sequences_for_training/IncN/fastas
mkdir -p plasmid_sequences_for_training/IncX1/fastas
mkdir -p output/amrfinder
mkdir -p output/integrated
mkdir -p output/figures

echo "  Directory structure ready."

# ── 3. AMRFinderPlus (optional) ────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  OPTIONAL: AMRFinderPlus Installation"
echo "============================================================"
echo ""

# Check if amrfinder is already available
if command -v amrfinder &>/dev/null; then
    echo "  AMRFinderPlus found: $(amrfinder --version 2>&1 | head -1)"
    echo "  Skipping installation."
elif command -v conda &>/dev/null || command -v mamba &>/dev/null; then
    echo "  conda/mamba detected."
    read -p "  Install AMRFinderPlus via conda? [y/N]: " install_amr
    if [[ "$install_amr" =~ ^[Yy]$ ]]; then
        if command -v mamba &>/dev/null; then
            mamba install -c bioconda -c conda-forge ncbi-amrfinderplus -y
        else
            conda install -c bioconda -c conda-forge ncbi-amrfinderplus -y
        fi
        amrfinder -u  # Update database
        echo "  AMRFinderPlus installed and database updated."
    else
        echo "  Skipping AMRFinderPlus installation."
        echo "  You can install it later with:"
        echo "    conda install -c bioconda -c conda-forge ncbi-amrfinderplus"
    fi
else
    echo "  AMRFinderPlus not found and conda/mamba not available."
    echo "  To install AMRFinderPlus, first install conda/mamba, then run:"
    echo "    conda install -c bioconda -c conda-forge ncbi-amrfinderplus"
    echo "    amrfinder -u"
fi

echo ""
echo "============================================================"
echo "  Setup Complete!"
echo "============================================================"
echo ""
echo "  Next steps:"
echo "  1. Place your FASTA files in:"
echo "       plasmid_sequences_for_training/IncFII/fastas/"
echo "       plasmid_sequences_for_training/IncN/fastas/"
echo "       plasmid_sequences_for_training/IncX1/fastas/"
echo ""
echo "  2. Run the full pipeline:"
echo "       bash run_all.sh"
echo ""
echo "  Or run individual steps:"
echo "       source .venv/bin/activate"
echo "       python assign_pLIN.py"
echo "       bash run_amrfinder_all.sh"
echo "       python integrate_pLIN_AMR.py"
echo "       python generate_figures.py"
echo ""
echo "============================================================"
