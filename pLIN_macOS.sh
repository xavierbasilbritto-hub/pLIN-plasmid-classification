#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════════
#  pLIN — All-in-One Setup & Launch for macOS
#  Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
#
#  Single file: installs everything + launches the pLIN web application.
#
#  Usage:
#    chmod +x pLIN_macOS.sh
#    ./pLIN_macOS.sh                  # Full install + launch
#    ./pLIN_macOS.sh --install        # Install only (no launch)
#    ./pLIN_macOS.sh --launch         # Launch only (skip install)
#    ./pLIN_macOS.sh --check          # Check environment only
#    ./pLIN_macOS.sh --uninstall      # Remove conda env + venv
#    ./pLIN_macOS.sh --docker         # Build & run via Docker
#
#  Works on: macOS Ventura 13+, Sonoma 14+, Sequoia 15+
#            Apple Silicon (M1/M2/M3/M4) and Intel Macs
# ═══════════════════════════════════════════════════════════════════════════════
set -euo pipefail

# ── Constants ─────────────────────────────────────────────────────────────────
APP_NAME="pLIN"
APP_VERSION="3.0.0"
ENV_NAME="pLIN_tools"
PYTHON_MIN="3.10"
STREAMLIT_PORT=8501
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_FILE="$SCRIPT_DIR/plin_app.py"
REQUIREMENTS="$SCRIPT_DIR/requirements.txt"

# Terminal colours
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'  # No colour

# ── Helper functions ──────────────────────────────────────────────────────────
banner() {
    echo ""
    echo -e "${CYAN}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${CYAN}║${NC}  ${BOLD}pLIN: Plasmid Lineage Identification Number System${NC}        ${CYAN}║${NC}"
    echo -e "${CYAN}║${NC}  Version ${APP_VERSION} — All-in-One Installer for macOS          ${CYAN}║${NC}"
    echo -e "${CYAN}║${NC}  Hierarchical Plasmid Classification + AMR Surveillance    ${CYAN}║${NC}"
    echo -e "${CYAN}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""
}

ok()   { echo -e "  ${GREEN}[OK]${NC}   $1"; }
warn() { echo -e "  ${YELLOW}[WARN]${NC} $1"; }
fail() { echo -e "  ${RED}[FAIL]${NC} $1"; }
info() { echo -e "  ${CYAN}[INFO]${NC} $1"; }
step() { echo -e "\n${BOLD}[$1/$2] $3${NC}"; }

version_ge() {
    # Returns 0 (true) if $1 >= $2 (version comparison)
    printf '%s\n%s' "$2" "$1" | sort -V -C
}

# ── Detect platform ──────────────────────────────────────────────────────────
detect_platform() {
    ARCH=$(uname -m)
    if [ "$ARCH" = "arm64" ]; then
        PLATFORM="macOS (Apple Silicon)"
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
    else
        PLATFORM="macOS (Intel)"
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    fi
    MAC_VERSION=$(sw_vers -productVersion 2>/dev/null || echo "Unknown")
    echo -e "  Platform:     ${BOLD}$PLATFORM${NC}"
    echo -e "  macOS:        $MAC_VERSION"
    echo -e "  Architecture: $ARCH"
}

# ── Find or install conda ────────────────────────────────────────────────────
find_conda() {
    # Check PATH first
    if command -v conda &>/dev/null; then
        CONDA_BIN=$(command -v conda)
        return 0
    fi
    # Search common locations
    for loc in \
        "$HOME/miniconda3/bin/conda" \
        "$HOME/miniforge3/bin/conda" \
        "$HOME/anaconda3/bin/conda" \
        "$HOME/mambaforge/bin/conda" \
        "$HOME/opt/miniconda3/bin/conda" \
        "/opt/homebrew/Caskroom/miniconda/base/bin/conda" \
        "/usr/local/Caskroom/miniconda/base/bin/conda"; do
        if [ -x "$loc" ]; then
            CONDA_BIN="$loc"
            return 0
        fi
    done
    return 1
}

install_conda() {
    info "Conda not found. Installing Miniconda..."
    echo ""
    INSTALLER="/tmp/miniconda_installer.sh"
    curl -fsSL "$MINICONDA_URL" -o "$INSTALLER"
    bash "$INSTALLER" -b -p "$HOME/miniconda3"
    rm -f "$INSTALLER"
    CONDA_BIN="$HOME/miniconda3/bin/conda"
    # Initialize for this shell
    eval "$($CONDA_BIN shell.bash hook)"
    ok "Miniconda installed to $HOME/miniconda3"
}

# ── Find Python 3.10+ ────────────────────────────────────────────────────────
find_python() {
    for cmd in python3.12 python3.11 python3.10 python3 python; do
        if command -v "$cmd" &>/dev/null; then
            local ver=$("$cmd" --version 2>&1 | awk '{print $2}')
            if version_ge "$ver" "$PYTHON_MIN"; then
                PYTHON_BIN=$(command -v "$cmd")
                PYTHON_VER="$ver"
                return 0
            fi
        fi
    done
    return 1
}

# ── Check required data files ────────────────────────────────────────────────
check_files() {
    local all_ok=true
    local files=(
        "$APP_FILE"
        "$REQUIREMENTS"
        "$SCRIPT_DIR/data/inc_classifier.npz"
        "$SCRIPT_DIR/data/inc_centroids.npz"
        "$SCRIPT_DIR/output/pLIN_assignments.tsv"
    )
    for f in "${files[@]}"; do
        if [ -f "$f" ]; then
            ok "Found $(basename "$f")"
        else
            fail "Missing: $f"
            all_ok=false
        fi
    done
    $all_ok
}

# ── Check Python packages ────────────────────────────────────────────────────
check_packages() {
    local missing=()
    local packages=(streamlit numpy pandas scipy Bio sklearn matplotlib seaborn plotly pptx requests joblib)
    local names=(streamlit numpy pandas scipy biopython scikit-learn matplotlib seaborn plotly python-pptx requests joblib)

    for i in "${!packages[@]}"; do
        if "$PYTHON_BIN" -c "import ${packages[$i]}" 2>/dev/null; then
            ok "Python package: ${names[$i]}"
        else
            warn "Missing: ${names[$i]}"
            missing+=("${names[$i]}")
        fi
    done
    MISSING_PACKAGES=("${missing[@]+"${missing[@]}"}")
}

# ── Check bioinformatics tools ────────────────────────────────────────────────
check_biotools() {
    local tools=(amrfinder mash fastANI minimap2 minced blastn prodigal mob_typer mlst)
    local descs=(
        "AMR gene detection"
        "Fast genome distance"
        "Average nucleotide identity"
        "Sequence alignment / SNP typing"
        "CRISPR spacer extraction"
        "CRISPR host prediction"
        "Gene prediction"
        "Plasmid mobility typing"
        "Multi-locus sequence typing"
    )
    BIOTOOLS_FOUND=0
    BIOTOOLS_TOTAL=${#tools[@]}

    for i in "${!tools[@]}"; do
        if command -v "${tools[$i]}" &>/dev/null; then
            ok "${tools[$i]} — ${descs[$i]}"
            ((BIOTOOLS_FOUND++))
        else
            # Also check conda envs
            local found=false
            for base in "$HOME/miniconda3" "$HOME/miniforge3" "$HOME/anaconda3" "$HOME/mambaforge"; do
                if [ -x "$base/envs/$ENV_NAME/bin/${tools[$i]}" ]; then
                    ok "${tools[$i]} — ${descs[$i]} (in conda env)"
                    ((BIOTOOLS_FOUND++))
                    found=true
                    break
                fi
            done
            if ! $found; then
                warn "${tools[$i]} — ${descs[$i]} (not found, optional)"
            fi
        fi
    done
}

# ── Install Python packages ──────────────────────────────────────────────────
install_python_packages() {
    step "$1" "$2" "Installing Python dependencies..."

    # Use conda env python if available
    if [ -n "${CONDA_ENV_PYTHON:-}" ] && [ -x "$CONDA_ENV_PYTHON" ]; then
        PYTHON_BIN="$CONDA_ENV_PYTHON"
    fi

    "$PYTHON_BIN" -m pip install --upgrade pip --quiet 2>/dev/null || true
    "$PYTHON_BIN" -m pip install --quiet -r "$REQUIREMENTS"
    ok "All Python packages installed"
}

# ── Install bioinformatics tools via conda ────────────────────────────────────
install_biotools() {
    step "$1" "$2" "Installing bioinformatics tools via conda..."

    if ! find_conda; then
        install_conda
    fi
    ok "Conda: $CONDA_BIN"

    # Initialize conda for this shell
    eval "$($CONDA_BIN shell.bash hook)" 2>/dev/null || true

    # Create or activate conda environment
    if "$CONDA_BIN" env list 2>/dev/null | grep -q "^${ENV_NAME} "; then
        info "Conda environment '$ENV_NAME' already exists"
    else
        info "Creating conda environment '$ENV_NAME' with Python 3.11..."
        "$CONDA_BIN" create -n "$ENV_NAME" python=3.11 -y --quiet 2>/dev/null
        ok "Conda environment created"
    fi

    conda activate "$ENV_NAME" 2>/dev/null || \
        eval "$($CONDA_BIN shell.bash hook)" && conda activate "$ENV_NAME"

    # Save env python path
    CONDA_ENV_PYTHON="$CONDA_PREFIX/bin/python3"

    # Configure channels
    conda config --add channels defaults 2>/dev/null || true
    conda config --add channels bioconda 2>/dev/null || true
    conda config --add channels conda-forge 2>/dev/null || true
    conda config --set channel_priority strict 2>/dev/null || true

    # Install bioinformatics tools
    local tools=(ncbi-amrfinderplus mash fastani minimap2 minced blast prodigal mlst)
    local installed=0
    local total=${#tools[@]}

    for tool in "${tools[@]}"; do
        info "Installing $tool..."
        if conda install -y -c bioconda -c conda-forge "$tool" --quiet 2>/dev/null; then
            ok "Installed $tool"
            ((installed++))
        else
            warn "Could not install $tool (optional)"
        fi
    done

    # MOBsuite via pip (not in conda)
    info "Installing MOBsuite..."
    if pip install mob_suite --quiet 2>/dev/null; then
        ok "Installed MOBsuite"
        ((installed++))
    else
        warn "Could not install MOBsuite (optional)"
    fi

    # Update AMRFinderPlus database
    if command -v amrfinder &>/dev/null; then
        info "Updating AMRFinderPlus database..."
        amrfinder --update 2>/dev/null || warn "AMRFinderPlus database update failed (run manually)"
    fi

    ok "Installed $installed/$((total+1)) bioinformatics tools"
}

# ── Create directory structure ────────────────────────────────────────────────
setup_directories() {
    mkdir -p "$SCRIPT_DIR/output/amrfinder"
    mkdir -p "$SCRIPT_DIR/output/crispr_evaluation"
    mkdir -p "$SCRIPT_DIR/output/figures"
    mkdir -p "$SCRIPT_DIR/output/manuscripts/docx"
}

# ── Docker mode ───────────────────────────────────────────────────────────────
docker_mode() {
    banner
    step 1 3 "Checking Docker..."

    if ! command -v docker &>/dev/null; then
        fail "Docker not found"
        info "Install Docker Desktop from: https://www.docker.com/products/docker-desktop/"
        exit 1
    fi
    ok "Docker found"

    step 2 3 "Building pLIN Docker image..."
    docker build -t plin:latest "$SCRIPT_DIR"
    ok "Docker image built: plin:latest"

    step 3 3 "Starting pLIN container..."
    docker stop plin-app 2>/dev/null || true
    docker rm plin-app 2>/dev/null || true
    docker run -d --name plin-app \
        -p ${STREAMLIT_PORT}:${STREAMLIT_PORT} \
        -v "$SCRIPT_DIR/data:/app/data" \
        -v "$SCRIPT_DIR/output:/app/output" \
        plin:latest
    ok "pLIN running in Docker"
    echo ""
    echo -e "  ${BOLD}Open in browser: http://localhost:${STREAMLIT_PORT}${NC}"
    echo ""
    info "Stop with:  docker stop plin-app"
    info "Logs:       docker logs -f plin-app"
}

# ── Uninstall ─────────────────────────────────────────────────────────────────
uninstall() {
    banner
    echo -e "  ${YELLOW}This will remove the pLIN conda environment and virtual environment.${NC}"
    echo -e "  ${YELLOW}Your data files and sequences will NOT be deleted.${NC}"
    echo ""
    read -p "  Continue? [y/N]: " confirm
    if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
        echo "  Cancelled."
        exit 0
    fi

    if find_conda; then
        eval "$($CONDA_BIN shell.bash hook)" 2>/dev/null || true
        if "$CONDA_BIN" env list 2>/dev/null | grep -q "^${ENV_NAME} "; then
            info "Removing conda environment '$ENV_NAME'..."
            "$CONDA_BIN" env remove -n "$ENV_NAME" -y 2>/dev/null || true
            ok "Conda environment removed"
        fi
    fi

    if [ -d "$SCRIPT_DIR/.venv" ]; then
        info "Removing virtual environment..."
        rm -rf "$SCRIPT_DIR/.venv"
        ok "Virtual environment removed"
    fi

    ok "Uninstall complete. Data files are preserved."
}

# ── Launch app ────────────────────────────────────────────────────────────────
launch_app() {
    step "$1" "$2" "Launching pLIN..."

    if [ ! -f "$APP_FILE" ]; then
        fail "Application file not found: $APP_FILE"
        exit 1
    fi

    # Try to activate conda env if it exists
    if find_conda; then
        eval "$($CONDA_BIN shell.bash hook)" 2>/dev/null || true
        if "$CONDA_BIN" env list 2>/dev/null | grep -q "^${ENV_NAME} "; then
            conda activate "$ENV_NAME" 2>/dev/null || true
        fi
    fi

    # Find streamlit
    local STREAMLIT=""
    if command -v streamlit &>/dev/null; then
        STREAMLIT="streamlit"
    elif [ -x "$CONDA_PREFIX/bin/streamlit" 2>/dev/null ]; then
        STREAMLIT="$CONDA_PREFIX/bin/streamlit"
    elif [ -x "$SCRIPT_DIR/.venv/bin/streamlit" ]; then
        STREAMLIT="$SCRIPT_DIR/.venv/bin/streamlit"
    fi

    if [ -z "$STREAMLIT" ]; then
        fail "Streamlit not found. Run: ./pLIN_macOS.sh --install"
        exit 1
    fi

    echo ""
    echo -e "  ${BOLD}╔══════════════════════════════════════════════════╗${NC}"
    echo -e "  ${BOLD}║  pLIN is starting at: http://localhost:${STREAMLIT_PORT}    ║${NC}"
    echo -e "  ${BOLD}║  Press Ctrl+C to stop the server.               ║${NC}"
    echo -e "  ${BOLD}╚══════════════════════════════════════════════════╝${NC}"
    echo ""

    # Open browser automatically on macOS
    ( sleep 3 && open "http://localhost:${STREAMLIT_PORT}" 2>/dev/null ) &

    "$STREAMLIT" run "$APP_FILE" \
        --server.headless true \
        --server.port "$STREAMLIT_PORT" \
        --server.address 0.0.0.0 \
        --browser.gatherUsageStats false \
        --server.maxUploadSize 200 \
        || true
}

# ── Check-only mode ───────────────────────────────────────────────────────────
check_mode() {
    banner
    detect_platform
    echo ""

    step 1 4 "Checking Python..."
    if find_python; then
        ok "Python $PYTHON_VER ($PYTHON_BIN)"
    else
        fail "Python $PYTHON_MIN+ not found"
    fi

    step 2 4 "Checking required files..."
    check_files || true

    step 3 4 "Checking Python packages..."
    if find_python; then
        check_packages
    else
        warn "Cannot check packages — Python not found"
    fi

    step 4 4 "Checking bioinformatics tools..."
    check_biotools

    # Summary
    echo ""
    echo -e "${BOLD}══════════════════════════════════════════════════════════════${NC}"
    echo -e "  ${BOLD}Environment Summary${NC}"
    echo -e "${BOLD}══════════════════════════════════════════════════════════════${NC}"
    if [ ${#MISSING_PACKAGES[@]} -eq 0 ] 2>/dev/null; then
        ok "Core Python packages: all installed"
    else
        warn "Missing ${#MISSING_PACKAGES[@]} Python packages"
    fi
    ok "Bioinformatics tools: $BIOTOOLS_FOUND/$BIOTOOLS_TOTAL available"
    echo ""
}

# ── Full install ──────────────────────────────────────────────────────────────
full_install() {
    local total_steps=5

    step 1 $total_steps "Checking environment..."
    detect_platform
    echo ""

    if find_python; then
        ok "Python $PYTHON_VER"
    else
        warn "Python $PYTHON_MIN+ not found in PATH"
        info "Will install via conda"
    fi

    check_files || {
        fail "Required data files are missing."
        fail "Ensure you have the complete pLIN distribution."
        exit 1
    }

    step 2 $total_steps "Setting up conda environment + bioinformatics tools..."
    install_biotools 2 $total_steps

    step 3 $total_steps "Installing Python packages in conda environment..."
    install_python_packages 3 $total_steps

    step 4 $total_steps "Setting up directories..."
    setup_directories
    ok "Output directories ready"

    ok "Installation complete!"
}

# ── Main ──────────────────────────────────────────────────────────────────────
main() {
    cd "$SCRIPT_DIR"

    # Parse arguments
    case "${1:-}" in
        --install)
            banner
            full_install
            echo ""
            echo -e "  ${BOLD}Installation complete.${NC}"
            echo -e "  Launch pLIN with: ${CYAN}./pLIN_macOS.sh --launch${NC}"
            echo -e "  Or directly:      ${CYAN}conda activate $ENV_NAME && streamlit run plin_app.py${NC}"
            echo ""
            ;;
        --launch)
            banner
            launch_app 1 1
            ;;
        --check)
            check_mode
            ;;
        --docker)
            docker_mode
            ;;
        --uninstall)
            uninstall
            ;;
        --help|-h)
            banner
            echo "  Usage: ./pLIN_macOS.sh [OPTION]"
            echo ""
            echo "  Options:"
            echo "    (no option)    Full install + launch (default)"
            echo "    --install      Install only (no launch)"
            echo "    --launch       Launch only (skip install)"
            echo "    --check        Check environment status"
            echo "    --docker       Build and run via Docker"
            echo "    --uninstall    Remove conda env + venv"
            echo "    --help         Show this help message"
            echo ""
            ;;
        *)
            banner
            full_install
            launch_app 5 5
            ;;
    esac
}

main "$@"
