#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════
# pLIN Tool — Installation Script (macOS / Linux)
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# ═══════════════════════════════════════════════════════════════════════════
set -e

echo "═══════════════════════════════════════════════════════════════"
echo " pLIN Tool — Installer for macOS / Linux"
echo "═══════════════════════════════════════════════════════════════"
echo ""

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ENV_NAME="pLIN_tools"

# ── Check for conda ──────────────────────────────────────────────────────
if command -v conda &>/dev/null; then
    echo "[OK] Conda found: $(conda --version)"
else
    echo "[!] Conda not found."
    echo "    Install Miniconda from: https://docs.conda.io/en/latest/miniconda.html"
    echo ""
    echo "    macOS:  brew install --cask miniconda"
    echo "    Linux:  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "            bash Miniconda3-latest-Linux-x86_64.sh"
    echo ""
    exit 1
fi

# ── Create / activate conda environment ──────────────────────────────────
echo ""
echo "[1/5] Setting up conda environment: $ENV_NAME"
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "      Environment '$ENV_NAME' already exists. Activating..."
else
    echo "      Creating new environment with Python 3.11..."
    conda create -n "$ENV_NAME" python=3.11 -y
fi

# Activate
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"
echo "      Python: $(python3 --version)"

# ── Install Python dependencies ──────────────────────────────────────────
echo ""
echo "[2/5] Installing Python dependencies..."
pip install --upgrade pip
pip install \
    streamlit \
    numpy \
    pandas \
    scipy \
    biopython \
    scikit-learn \
    matplotlib \
    seaborn \
    plotly \
    python-pptx \
    requests

echo "      Core Python packages installed."

# ── Install bioinformatics tools via conda ───────────────────────────────
echo ""
echo "[3/5] Installing bioinformatics tools (conda-forge + bioconda)..."

# Ensure bioconda channel is available
conda config --add channels bioconda 2>/dev/null || true
conda config --add channels conda-forge 2>/dev/null || true

echo "      Installing AMRFinderPlus..."
conda install -c bioconda -c conda-forge ncbi-amrfinderplus -y 2>/dev/null && \
    echo "      [OK] AMRFinderPlus installed" || \
    echo "      [SKIP] AMRFinderPlus install failed (optional)"

echo "      Installing Mash..."
conda install -c bioconda mash -y 2>/dev/null && \
    echo "      [OK] Mash installed" || \
    echo "      [SKIP] Mash install failed (optional)"

echo "      Installing FastANI..."
conda install -c bioconda fastani -y 2>/dev/null && \
    echo "      [OK] FastANI installed" || \
    echo "      [SKIP] FastANI install failed (optional)"

echo "      Installing minimap2..."
conda install -c bioconda minimap2 -y 2>/dev/null && \
    echo "      [OK] minimap2 installed" || \
    echo "      [SKIP] minimap2 install failed (optional)"

echo "      Installing MinCED..."
conda install -c bioconda minced -y 2>/dev/null && \
    echo "      [OK] MinCED installed" || \
    echo "      [SKIP] MinCED install failed (optional)"

echo "      Installing BLAST+..."
conda install -c bioconda blast -y 2>/dev/null && \
    echo "      [OK] BLAST+ installed" || \
    echo "      [SKIP] BLAST+ install failed (optional)"

echo "      Installing Prodigal..."
conda install -c bioconda prodigal -y 2>/dev/null && \
    echo "      [OK] Prodigal installed" || \
    echo "      [SKIP] Prodigal install failed (optional)"

# ── Install MOBsuite via pip ─────────────────────────────────────────────
echo ""
echo "[4/5] Installing MOBsuite..."
pip install mob_suite 2>/dev/null && \
    echo "      [OK] MOBsuite installed" || \
    echo "      [SKIP] MOBsuite install failed (optional)"

# ── Update AMRFinderPlus database ────────────────────────────────────────
echo ""
echo "[5/5] Updating AMRFinderPlus database..."
amrfinder --update 2>/dev/null && \
    echo "      [OK] AMRFinderPlus database updated" || \
    echo "      [SKIP] AMRFinderPlus database update failed (run manually: amrfinder --update)"

# ── Summary ──────────────────────────────────────────────────────────────
echo ""
echo "═══════════════════════════════════════════════════════════════"
echo " Installation Complete!"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo " To launch pLIN:"
echo "   conda activate $ENV_NAME"
echo "   streamlit run $SCRIPT_DIR/plin_app.py"
echo ""
echo " The app will open at: http://localhost:8501"
echo ""
echo " Optional: Install Ollama for AI chatbot (Bacterial Buddy):"
echo "   macOS:  brew install ollama && ollama pull llama3.2"
echo "   Linux:  curl -fsSL https://ollama.com/install.sh | sh && ollama pull llama3.2"
echo ""
echo "═══════════════════════════════════════════════════════════════"
