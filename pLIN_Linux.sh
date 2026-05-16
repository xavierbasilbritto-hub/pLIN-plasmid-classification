#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════════
#  pLIN — All-in-One Setup & Launch for Linux
#  Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
#
#  Single file: installs everything + launches the pLIN web application.
#
#  Usage:
#    chmod +x pLIN_Linux.sh
#    ./pLIN_Linux.sh                  # Full install + launch
#    ./pLIN_Linux.sh --install        # Install only (no launch)
#    ./pLIN_Linux.sh --launch         # Launch only (skip install)
#    ./pLIN_Linux.sh --check          # Check environment only
#    ./pLIN_Linux.sh --uninstall      # Remove conda env + venv
#    ./pLIN_Linux.sh --docker         # Build & run via Docker
#
#  Works on: Ubuntu 20.04+, Debian 11+, Fedora 36+, CentOS 8+,
#            Arch Linux, RHEL 8+, openSUSE 15+
#            x86_64 and aarch64 (ARM64)
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
NC='\033[0m'

# ── Helper functions ──────────────────────────────────────────────────────────
banner() {
    echo ""
    echo -e "${CYAN}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${CYAN}║${NC}  ${BOLD}pLIN: Plasmid Lineage Identification Number System${NC}        ${CYAN}║${NC}"
    echo -e "${CYAN}║${NC}  Version ${APP_VERSION} — All-in-One Installer for Linux           ${CYAN}║${NC}"
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
    printf '%s\n%s' "$2" "$1" | sort -V -C
}

# ── Detect platform ──────────────────────────────────────────────────────────
detect_platform() {
    ARCH=$(uname -m)
    DISTRO="Unknown Linux"
    DISTRO_ID=""
    PKG_MGR=""

    # Detect distribution
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        DISTRO="$NAME $VERSION_ID"
        DISTRO_ID="$ID"
    elif [ -f /etc/redhat-release ]; then
        DISTRO=$(cat /etc/redhat-release)
        DISTRO_ID="rhel"
    fi

    # Detect package manager
    if command -v apt-get &>/dev/null; then
        PKG_MGR="apt"
    elif command -v dnf &>/dev/null; then
        PKG_MGR="dnf"
    elif command -v yum &>/dev/null; then
        PKG_MGR="yum"
    elif command -v pacman &>/dev/null; then
        PKG_MGR="pacman"
    elif command -v zypper &>/dev/null; then
        PKG_MGR="zypper"
    fi

    # Miniconda URL
    if [ "$ARCH" = "aarch64" ] || [ "$ARCH" = "arm64" ]; then
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh"
    else
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    fi

    echo -e "  Distribution: ${BOLD}$DISTRO${NC}"
    echo -e "  Architecture: $ARCH"
    echo -e "  Pkg manager:  ${PKG_MGR:-none detected}"
    echo -e "  Kernel:       $(uname -r)"
}

# ── Install system dependencies ───────────────────────────────────────────────
install_system_deps() {
    info "Checking system dependencies..."

    # Check if curl/wget are available
    if ! command -v curl &>/dev/null && ! command -v wget &>/dev/null; then
        warn "Neither curl nor wget found. Installing..."
        case "$PKG_MGR" in
            apt)    sudo apt-get update -qq && sudo apt-get install -y -qq curl wget build-essential default-jre-headless ;;
            dnf)    sudo dnf install -y curl wget gcc gcc-c++ java-11-openjdk-headless ;;
            yum)    sudo yum install -y curl wget gcc gcc-c++ java-11-openjdk-headless ;;
            pacman) sudo pacman -Sy --noconfirm curl wget base-devel jre-openjdk-headless ;;
            zypper) sudo zypper install -y curl wget gcc gcc-c++ java-11-openjdk-headless ;;
            *)      fail "Cannot install system deps — unknown package manager"; exit 1 ;;
        esac
    fi

    # Java is needed for MinCED
    if ! command -v java &>/dev/null; then
        warn "Java not found (needed for MinCED). Installing..."
        case "$PKG_MGR" in
            apt)    sudo apt-get install -y -qq default-jre-headless 2>/dev/null || true ;;
            dnf)    sudo dnf install -y java-11-openjdk-headless 2>/dev/null || true ;;
            yum)    sudo yum install -y java-11-openjdk-headless 2>/dev/null || true ;;
            pacman) sudo pacman -Sy --noconfirm jre-openjdk-headless 2>/dev/null || true ;;
            zypper) sudo zypper install -y java-11-openjdk-headless 2>/dev/null || true ;;
            *)      warn "Install Java manually for MinCED support" ;;
        esac
    fi

    ok "System dependencies checked"
}

# ── Find or install conda ────────────────────────────────────────────────────
find_conda() {
    if command -v conda &>/dev/null; then
        CONDA_BIN=$(command -v conda)
        return 0
    fi
    for loc in \
        "$HOME/miniconda3/bin/conda" \
        "$HOME/miniforge3/bin/conda" \
        "$HOME/anaconda3/bin/conda" \
        "$HOME/mambaforge/bin/conda" \
        "/opt/conda/bin/conda"; do
        if [ -x "$loc" ]; then
            CONDA_BIN="$loc"
            return 0
        fi
    done
    return 1
}

install_conda() {
    info "Conda not found. Installing Miniconda..."
    local INSTALLER="/tmp/miniconda_installer.sh"

    if command -v curl &>/dev/null; then
        curl -fsSL "$MINICONDA_URL" -o "$INSTALLER"
    elif command -v wget &>/dev/null; then
        wget -q "$MINICONDA_URL" -O "$INSTALLER"
    else
        fail "Neither curl nor wget available. Install one and try again."
        exit 1
    fi

    bash "$INSTALLER" -b -p "$HOME/miniconda3"
    rm -f "$INSTALLER"
    CONDA_BIN="$HOME/miniconda3/bin/conda"
    eval "$($CONDA_BIN shell.bash hook)"
    ok "Miniconda installed to $HOME/miniconda3"
    info "Run 'conda init bash' later to make conda available in new terminals"
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

    eval "$($CONDA_BIN shell.bash hook)" 2>/dev/null || true

    if "$CONDA_BIN" env list 2>/dev/null | grep -q "^${ENV_NAME} "; then
        info "Conda environment '$ENV_NAME' already exists"
    else
        info "Creating conda environment '$ENV_NAME' with Python 3.11..."
        "$CONDA_BIN" create -n "$ENV_NAME" python=3.11 -y --quiet 2>/dev/null
        ok "Conda environment created"
    fi

    conda activate "$ENV_NAME" 2>/dev/null || \
        eval "$($CONDA_BIN shell.bash hook)" && conda activate "$ENV_NAME"

    CONDA_ENV_PYTHON="$CONDA_PREFIX/bin/python3"

    conda config --add channels defaults 2>/dev/null || true
    conda config --add channels bioconda 2>/dev/null || true
    conda config --add channels conda-forge 2>/dev/null || true
    conda config --set channel_priority strict 2>/dev/null || true

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

    # MOBsuite via pip
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
        amrfinder --update 2>/dev/null || warn "AMRFinderPlus database update failed"
    fi

    ok "Installed $installed/$((total+1)) bioinformatics tools"
}

# ── Create directory structure ────────────────────────────────────────────────
setup_directories() {
    mkdir -p "$SCRIPT_DIR/output/amrfinder"
    mkdir -p "$SCRIPT_DIR/output/crispr_evaluation"
    mkdir -p "$SCRIPT_DIR/output/figures"
    mkdir -p "$SCRIPT_DIR/output/manuscripts/docx"
    mkdir -p "$SCRIPT_DIR/output/mge_detection"
}

# ── IS reference database setup ──────────────────────────────────────────────
setup_is_database() {
    local IS_DIR="$SCRIPT_DIR/output/mge_detection"
    local IS_FASTA="$IS_DIR/is_reference_sequences.fasta"
    local IS_DB="$IS_DIR/is_reference_sequences.fasta.ndb"

    # Skip if database already exists
    if [ -f "$IS_DB" ]; then
        ok "IS reference BLAST database already exists"
        return 0
    fi

    info "Setting up IS element reference database (25 IS families)..."

    # Find python in conda env
    local PY=""
    if [ -n "${CONDA_ENV_PYTHON:-}" ] && [ -x "$CONDA_ENV_PYTHON" ]; then
        PY="$CONDA_ENV_PYTHON"
    elif [ -n "${PYTHON_BIN:-}" ]; then
        PY="$PYTHON_BIN"
    else
        warn "Python not available — skipping IS database setup"
        return 0
    fi

    # Download IS reference sequences from NCBI via BioPython
    "$PY" -c "
import sys, os
try:
    from Bio import Entrez, SeqIO
except ImportError:
    print('  [WARN] BioPython not available — skipping IS database download')
    sys.exit(0)

Entrez.email = 'plin_tool@example.com'
accessions = {
    'IS26': 'X00011.1', 'ISEcp1': 'AJ242809.1', 'IS1': 'J01730.1',
    'IS903': 'M17148.1', 'IS6100': 'M95400.1', 'ISKpn26': 'KF914891.1',
    'IS5': 'X02311.1', 'IS3': 'X02180.1', 'IS4321': 'AJ245418.1',
    'IS15': 'X01840.1', 'IS10': 'J01830.1', 'IS2': 'J01733.1',
    'IS4': 'V00029.1', 'IS30': 'X00792.1', 'IS66': 'X53365.1',
    'IS110': 'M21395.1', 'ISPa': 'AF261825.1', 'ISAba': 'AY758396.1',
    'IS256': 'M18086.1', 'IS257': 'U40412.1', 'IS16': 'AF053365.1',
    'ISEnfa': 'AF162694.1', 'IS1216': 'L40841.1', 'IS1251': 'X83579.1',
    'Tn916': 'U09422.1'
}

out_path = '$IS_FASTA'
os.makedirs(os.path.dirname(out_path), exist_ok=True)
downloaded = 0
with open(out_path, 'w') as fh:
    for name, acc in accessions.items():
        try:
            handle = Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text')
            record = SeqIO.read(handle, 'fasta')
            handle.close()
            record.id = name
            record.description = f'{name} ({acc})'
            SeqIO.write(record, fh, 'fasta')
            downloaded += 1
        except Exception as e:
            print(f'  [WARN] Could not download {name} ({acc}): {e}')

print(f'  [OK]   Downloaded {downloaded}/25 IS reference sequences')
" 2>/dev/null

    if [ $? -ne 0 ]; then
        warn "IS reference download failed (NCBI may be unreachable)"
        return 0
    fi

    # Build BLAST database
    if [ -f "$IS_FASTA" ] && command -v makeblastdb &>/dev/null; then
        makeblastdb -in "$IS_FASTA" -dbtype nucl -parse_seqids \
            -title "IS_reference_25families" -out "$IS_FASTA" &>/dev/null
        if [ $? -eq 0 ]; then
            ok "BLAST database built for IS references"
        else
            warn "makeblastdb failed — IS detection may not work"
        fi
    elif [ -f "$IS_FASTA" ]; then
        warn "makeblastdb not found — BLAST database not built"
        info "IS detection requires BLAST+ (installed via conda)"
    fi
}

# ── Verify BLAST for contig classification and IS detection ──────────────────
verify_blast() {
    if command -v blastn &>/dev/null; then
        local ver=$(blastn -version 2>&1 | head -1)
        ok "BLAST+ available: $ver"
    else
        # Check conda env
        local found=false
        for base in "$HOME/miniconda3" "$HOME/miniforge3" "$HOME/anaconda3" "$HOME/mambaforge"; do
            if [ -x "$base/envs/$ENV_NAME/bin/blastn" ]; then
                local ver=$("$base/envs/$ENV_NAME/bin/blastn" -version 2>&1 | head -1)
                ok "BLAST+ available (in conda env): $ver"
                found=true
                break
            fi
        done
        if ! $found; then
            warn "BLAST+ not found — IS element detection will not be available"
            info "Contig classification (plasmid vs chromosome) works without BLAST"
        fi
    fi
}

# ── Docker mode ───────────────────────────────────────────────────────────────
docker_mode() {
    banner
    step 1 3 "Checking Docker..."

    if ! command -v docker &>/dev/null; then
        fail "Docker not found"
        echo ""
        info "Install Docker:"
        info "  Ubuntu/Debian: sudo apt install docker.io docker-compose"
        info "  Fedora/RHEL:   sudo dnf install docker docker-compose"
        info "  Arch:          sudo pacman -S docker docker-compose"
        info "  Or:            https://docs.docker.com/engine/install/"
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

    # Try to activate conda env
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
    elif [ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/bin/streamlit" ]; then
        STREAMLIT="$CONDA_PREFIX/bin/streamlit"
    elif [ -x "$SCRIPT_DIR/.venv/bin/streamlit" ]; then
        STREAMLIT="$SCRIPT_DIR/.venv/bin/streamlit"
    fi

    if [ -z "$STREAMLIT" ]; then
        fail "Streamlit not found. Run: ./pLIN_Linux.sh --install"
        exit 1
    fi

    echo ""
    echo -e "  ${BOLD}╔══════════════════════════════════════════════════╗${NC}"
    echo -e "  ${BOLD}║  pLIN is starting at: http://localhost:${STREAMLIT_PORT}    ║${NC}"
    echo -e "  ${BOLD}║  Press Ctrl+C to stop the server.               ║${NC}"
    echo -e "  ${BOLD}╚══════════════════════════════════════════════════╝${NC}"
    echo ""

    # Open browser (try multiple methods for Linux)
    ( sleep 3 && {
        xdg-open "http://localhost:${STREAMLIT_PORT}" 2>/dev/null || \
        sensible-browser "http://localhost:${STREAMLIT_PORT}" 2>/dev/null || \
        x-www-browser "http://localhost:${STREAMLIT_PORT}" 2>/dev/null || \
        true
    } ) &

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

    step 1 5 "Checking Python..."
    if find_python; then
        ok "Python $PYTHON_VER ($PYTHON_BIN)"
    else
        fail "Python $PYTHON_MIN+ not found"
    fi

    step 2 5 "Checking required files..."
    check_files || true

    step 3 5 "Checking Python packages..."
    if find_python; then
        check_packages
    else
        warn "Cannot check packages — Python not found"
    fi

    step 4 5 "Checking bioinformatics tools..."
    check_biotools

    step 5 5 "Checking IS reference database + BLAST..."
    local IS_DB="$SCRIPT_DIR/output/mge_detection/is_reference_sequences.fasta.ndb"
    if [ -f "$IS_DB" ]; then
        ok "IS reference BLAST database found"
    else
        warn "IS reference database not built (run --install to set up)"
    fi
    verify_blast

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
    local total_steps=6

    step 1 $total_steps "Checking environment..."
    detect_platform
    echo ""
    install_system_deps

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

    step 5 $total_steps "Setting up IS element reference database + BLAST verification..."
    setup_is_database
    verify_blast

    ok "Installation complete!"
}

# ── Main ──────────────────────────────────────────────────────────────────────
main() {
    cd "$SCRIPT_DIR"

    case "${1:-}" in
        --install)
            banner
            full_install
            echo ""
            echo -e "  ${BOLD}Installation complete.${NC}"
            echo -e "  Launch pLIN with: ${CYAN}./pLIN_Linux.sh --launch${NC}"
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
            echo "  Usage: ./pLIN_Linux.sh [OPTION]"
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
            echo "  Supported distributions:"
            echo "    Ubuntu 20.04+, Debian 11+, Fedora 36+, CentOS 8+,"
            echo "    RHEL 8+, Arch Linux, openSUSE 15+ (x86_64 and ARM64)"
            echo ""
            ;;
        *)
            banner
            full_install
            launch_app 6 6
            ;;
    esac
}

main "$@"
