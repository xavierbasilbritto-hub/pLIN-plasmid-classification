#!/bin/bash
# ============================================================
#  pLIN Launcher for Linux
#  Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
#  Run: chmod +x launch_pLIN.sh && ./launch_pLIN.sh
# ============================================================
set -e

echo ""
echo " ============================================"
echo "  pLIN: Plasmid Life Identification Number"
echo "  Hierarchical Plasmid Classification System"
echo " ============================================"
echo ""

# Change to script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
echo "[OK] Working directory: $SCRIPT_DIR"

# Find Python 3.10+
PYTHON=""
for cmd in python3 python; do
    if command -v "$cmd" &>/dev/null; then
        version=$("$cmd" --version 2>&1 | awk '{print $2}')
        major=$(echo "$version" | cut -d. -f1)
        minor=$(echo "$version" | cut -d. -f2)
        if [ "$major" -ge 3 ] && [ "$minor" -ge 10 ]; then
            PYTHON="$cmd"
            echo "[OK] Found $cmd $version"
            break
        fi
    fi
done

if [ -z "$PYTHON" ]; then
    echo "[ERROR] Python 3.10+ is required but not found."
    echo ""
    echo "Install Python via your package manager:"
    echo "  Ubuntu/Debian: sudo apt install python3 python3-venv python3-pip"
    echo "  Fedora/RHEL:   sudo dnf install python3 python3-pip"
    echo "  Arch:          sudo pacman -S python python-pip"
    echo ""
    exit 1
fi

# Check for python3-venv (common issue on Debian/Ubuntu)
if ! "$PYTHON" -m venv --help &>/dev/null; then
    echo "[ERROR] python3-venv is not installed."
    echo "  Install via: sudo apt install python3-venv"
    exit 1
fi

# Create virtual environment if needed
if [ ! -f ".venv/bin/activate" ]; then
    echo ""
    echo "[SETUP] Creating virtual environment..."
    "$PYTHON" -m venv .venv
    echo "[OK] Virtual environment created."
fi

# Activate virtual environment
source .venv/bin/activate
echo "[OK] Virtual environment activated."

# Install dependencies
echo ""
echo "[SETUP] Installing dependencies (this may take a moment on first run)..."
pip install --quiet --upgrade pip
pip install --quiet -r requirements.txt
echo "[OK] All dependencies installed."

# Check for AMRFinderPlus
if command -v amrfinder &>/dev/null; then
    echo "[OK] AMRFinderPlus detected — AMR analysis available."
else
    echo "[INFO] AMRFinderPlus not found — pLIN will run without AMR analysis."
    echo "       Install via: conda install -c bioconda -c conda-forge ncbi-amrfinderplus"
fi

# Launch
echo ""
echo "============================================"
echo " Launching pLIN web application..."
echo " The app will open in your default browser."
echo " Press Ctrl+C in this window to stop."
echo "============================================"
echo ""
streamlit run plin_app.py --server.headless false --browser.gatherUsageStats false
