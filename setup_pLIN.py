#!/usr/bin/env python3
"""
pLIN Tool — Cross-Platform Setup & Launch Script
=================================================
Single file that installs dependencies and launches pLIN on any OS.

Usage:
    python3 setup_pLIN.py              # Install + launch
    python3 setup_pLIN.py --install    # Install only (no launch)
    python3 setup_pLIN.py --launch     # Launch only (skip install)
    python3 setup_pLIN.py --check      # Check environment only
    python3 setup_pLIN.py --docker     # Build and run Docker container

Works on: macOS (Intel/Apple Silicon), Linux (Ubuntu/Debian/Fedora/Arch), Windows 10/11

Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
"""

import os
import sys
import platform
import subprocess
import shutil
import argparse
import json
from pathlib import Path

# ── Constants ────────────────────────────────────────────────────────────
APP_NAME = "pLIN"
APP_VERSION = "2.1.0"
ENV_NAME = "pLIN_tools"
MIN_PYTHON = (3, 10)
STREAMLIT_PORT = 8501

SCRIPT_DIR = Path(__file__).resolve().parent
APP_FILE = SCRIPT_DIR / "plin_app.py"
REQUIREMENTS = SCRIPT_DIR / "requirements.txt"
DATA_DIR = SCRIPT_DIR / "data"
OUTPUT_DIR = SCRIPT_DIR / "output"

# Core data files required for pLIN to function
REQUIRED_FILES = [
    APP_FILE,
    REQUIREMENTS,
    DATA_DIR / "inc_classifier.npz",
    DATA_DIR / "inc_centroids.npz",
    OUTPUT_DIR / "pLIN_assignments.tsv",
]

# Python packages (core)
CORE_PACKAGES = [
    "streamlit", "numpy", "pandas", "scipy", "biopython",
    "scikit-learn", "matplotlib", "seaborn", "plotly",
    "python-pptx", "requests", "joblib",
]

# Bioinformatics tools (optional, installed via conda)
BIOTOOLS = {
    "amrfinder":    {"conda": "ncbi-amrfinderplus", "desc": "AMR gene detection"},
    "mash":         {"conda": "mash",               "desc": "Fast genome distance estimation"},
    "fastANI":      {"conda": "fastani",             "desc": "Average nucleotide identity"},
    "minimap2":     {"conda": "minimap2",            "desc": "Sequence alignment / SNP typing"},
    "minced":       {"conda": "minced",              "desc": "CRISPR spacer extraction"},
    "blastn":       {"conda": "blast",               "desc": "CRISPR host prediction"},
    "makeblastdb":  {"conda": "blast",               "desc": "BLAST database creation"},
    "prodigal":     {"conda": "prodigal",            "desc": "Gene prediction"},
}

# ── Colours for terminal output ──────────────────────────────────────────
class C:
    """ANSI colour codes (disabled on Windows CMD without VT support)."""
    _enabled = sys.stdout.isatty() and os.name != "nt"
    BOLD  = "\033[1m"   if _enabled else ""
    GREEN = "\033[92m"  if _enabled else ""
    YELLOW= "\033[93m"  if _enabled else ""
    RED   = "\033[91m"  if _enabled else ""
    CYAN  = "\033[96m"  if _enabled else ""
    RESET = "\033[0m"   if _enabled else ""


def banner():
    print(f"""
{C.CYAN}{'='*60}
  {C.BOLD}pLIN: Plasmid Life Identification Number System{C.RESET}{C.CYAN}
  Version {APP_VERSION} — Cross-Platform Setup
  Hierarchical Plasmid Classification + AMR Surveillance
{'='*60}{C.RESET}
""")


def ok(msg):
    print(f"  {C.GREEN}[OK]{C.RESET}   {msg}")

def warn(msg):
    print(f"  {C.YELLOW}[WARN]{C.RESET} {msg}")

def fail(msg):
    print(f"  {C.RED}[FAIL]{C.RESET} {msg}")

def info(msg):
    print(f"  {C.CYAN}[INFO]{C.RESET} {msg}")

def step(n, total, msg):
    print(f"\n{C.BOLD}[{n}/{total}] {msg}{C.RESET}")


# ── Platform detection ───────────────────────────────────────────────────
def detect_platform():
    """Detect OS and architecture."""
    system = platform.system().lower()
    machine = platform.machine().lower()
    is_arm = machine in ("arm64", "aarch64")

    if system == "darwin":
        os_name = "macOS (Apple Silicon)" if is_arm else "macOS (Intel)"
    elif system == "linux":
        os_name = f"Linux ({machine})"
    elif system == "windows":
        os_name = "Windows"
    else:
        os_name = f"{system} ({machine})"

    return {"system": system, "machine": machine, "is_arm": is_arm, "os_name": os_name}


# ── Command helpers ──────────────────────────────────────────────────────
def run_cmd(cmd, check=True, capture=True, timeout=600):
    """Run a shell command and return (success, stdout)."""
    try:
        r = subprocess.run(
            cmd, shell=isinstance(cmd, str), check=check,
            capture_output=capture, text=True, timeout=timeout,
        )
        return True, r.stdout.strip() if capture else ""
    except subprocess.CalledProcessError as e:
        return False, e.stderr.strip() if capture else str(e)
    except subprocess.TimeoutExpired:
        return False, "Command timed out"
    except FileNotFoundError:
        return False, "Command not found"


def which(tool):
    """Check if a tool is on PATH."""
    return shutil.which(tool)


def find_conda():
    """Find conda binary across common locations."""
    if which("conda"):
        return which("conda")

    # Search common locations
    home = Path.home()
    candidates = [
        home / "miniconda3" / "bin" / "conda",
        home / "miniforge3" / "bin" / "conda",
        home / "anaconda3" / "bin" / "conda",
        home / "mambaforge" / "bin" / "conda",
        Path("/opt/conda/bin/conda"),
    ]
    if platform.system() == "Windows":
        candidates += [
            home / "miniconda3" / "Scripts" / "conda.exe",
            home / "Miniconda3" / "Scripts" / "conda.exe",
            home / "Anaconda3" / "Scripts" / "conda.exe",
        ]
    for c in candidates:
        if c.exists():
            return str(c)
    return None


# ── Check environment ────────────────────────────────────────────────────
def check_python():
    """Verify Python version."""
    v = sys.version_info
    if v >= MIN_PYTHON:
        ok(f"Python {v.major}.{v.minor}.{v.micro}")
        return True
    else:
        fail(f"Python {v.major}.{v.minor}.{v.micro} — need {MIN_PYTHON[0]}.{MIN_PYTHON[1]}+")
        return False


def check_files():
    """Verify required data files exist."""
    all_ok = True
    for f in REQUIRED_FILES:
        if f.exists():
            ok(f"Found {f.name}")
        else:
            fail(f"Missing: {f}")
            all_ok = False
    return all_ok


def check_packages():
    """Verify Python packages are importable."""
    missing = []
    for pkg in CORE_PACKAGES:
        import_name = pkg.replace("-", "_")
        if import_name == "scikit_learn":
            import_name = "sklearn"
        if import_name == "python_pptx":
            import_name = "pptx"
        if import_name == "biopython":
            import_name = "Bio"
        try:
            __import__(import_name)
            ok(f"Python package: {pkg}")
        except ImportError:
            warn(f"Missing Python package: {pkg}")
            missing.append(pkg)
    return missing


def check_biotools():
    """Check which bioinformatics tools are available."""
    found = {}
    missing = {}
    for tool, meta in BIOTOOLS.items():
        path = which(tool)
        if path:
            found[tool] = path
            ok(f"{tool} — {meta['desc']}")
        else:
            missing[tool] = meta
            warn(f"{tool} not found — {meta['desc']} (optional)")
    return found, missing


# ── Installation ─────────────────────────────────────────────────────────
def install_python_packages():
    """Install Python dependencies via pip."""
    step(2, 4, "Installing Python dependencies...")
    success, output = run_cmd(
        [sys.executable, "-m", "pip", "install", "--upgrade", "pip"],
        check=False, capture=True,
    )
    success, output = run_cmd(
        [sys.executable, "-m", "pip", "install", "-r", str(REQUIREMENTS)],
        check=False, capture=True, timeout=300,
    )
    if success:
        ok("All Python packages installed")
        return True
    else:
        fail(f"pip install failed: {output[:200]}")
        info("Try manually: pip install -r requirements.txt")
        return False


def install_biotools_conda():
    """Install bioinformatics tools via conda (optional)."""
    step(3, 4, "Installing bioinformatics tools (optional)...")

    conda = find_conda()
    if not conda:
        warn("Conda not found — skipping bioinformatics tool installation")
        info("Core pLIN features (classification, pLIN assignment) work without these tools")
        info("Optional features requiring external tools: AMR detection, CRISPR analysis, ANI validation")
        info("Install conda from: https://docs.conda.io/en/latest/miniconda.html")
        return False

    ok(f"Conda found: {conda}")

    # Configure channels
    for ch in ["conda-forge", "bioconda", "defaults"]:
        run_cmd([conda, "config", "--add", "channels", ch], check=False)

    # Deduplicate conda packages (blast covers both blastn and makeblastdb)
    unique_pkgs = set()
    for meta in BIOTOOLS.values():
        unique_pkgs.add(meta["conda"])

    installed = 0
    for pkg in sorted(unique_pkgs):
        info(f"Installing {pkg}...")
        success, output = run_cmd(
            [conda, "install", "-y", "-c", "bioconda", "-c", "conda-forge", pkg],
            check=False, capture=True, timeout=300,
        )
        if success:
            ok(f"Installed {pkg}")
            installed += 1
        else:
            warn(f"Could not install {pkg} (optional)")

    # Update AMRFinderPlus database
    if which("amrfinder"):
        info("Updating AMRFinderPlus database...")
        run_cmd(["amrfinder", "--update"], check=False, capture=True, timeout=120)

    ok(f"Installed {installed}/{len(unique_pkgs)} bioinformatics tools")
    return True


# ── Docker ───────────────────────────────────────────────────────────────
def docker_build_and_run():
    """Build and run pLIN via Docker."""
    banner()
    step(1, 3, "Checking Docker...")

    if not which("docker"):
        fail("Docker not found")
        info("Install Docker Desktop from: https://www.docker.com/products/docker-desktop/")
        sys.exit(1)

    ok("Docker found")

    step(2, 3, "Building pLIN Docker image (this may take 5-10 minutes on first run)...")
    success, output = run_cmd(
        "docker build -t plin:latest .",
        check=False, capture=False, timeout=600,
    )
    if not success:
        fail("Docker build failed")
        sys.exit(1)
    ok("Docker image built: plin:latest")

    step(3, 3, "Starting pLIN container...")
    # Stop existing container if running
    run_cmd("docker stop plin-app 2>/dev/null", check=False)
    run_cmd("docker rm plin-app 2>/dev/null", check=False)

    success, output = run_cmd(
        f"docker run -d --name plin-app -p {STREAMLIT_PORT}:{STREAMLIT_PORT} "
        f"-v {SCRIPT_DIR}/data:/app/data "
        f"-v {SCRIPT_DIR}/output:/app/output "
        "plin:latest",
        check=False, capture=True,
    )
    if success:
        ok(f"pLIN running in Docker container")
        print(f"\n{C.BOLD}  Open in browser: http://localhost:{STREAMLIT_PORT}{C.RESET}\n")
        info("Stop with:  docker stop plin-app")
        info("Logs:       docker logs -f plin-app")
    else:
        fail(f"Failed to start container: {output[:200]}")


# ── Launch ───────────────────────────────────────────────────────────────
def launch():
    """Launch the Streamlit app."""
    step(4, 4, "Launching pLIN...")

    if not APP_FILE.exists():
        fail(f"Application file not found: {APP_FILE}")
        sys.exit(1)

    print(f"\n{C.BOLD}  pLIN is starting at: http://localhost:{STREAMLIT_PORT}{C.RESET}")
    print(f"  Press Ctrl+C to stop.\n")

    try:
        subprocess.run(
            [
                sys.executable, "-m", "streamlit", "run", str(APP_FILE),
                "--server.headless", "true",
                f"--server.port={STREAMLIT_PORT}",
                "--server.address=0.0.0.0",
                "--browser.gatherUsageStats=false",
                "--server.maxUploadSize=200",
            ],
            check=True,
        )
    except KeyboardInterrupt:
        print(f"\n{C.CYAN}pLIN stopped.{C.RESET}")
    except subprocess.CalledProcessError:
        fail("Streamlit failed to start")
        info("Try running directly: streamlit run plin_app.py")


# ── Environment report ───────────────────────────────────────────────────
def environment_report():
    """Generate a full environment report."""
    banner()
    plat = detect_platform()
    print(f"  Platform: {plat['os_name']}")
    print(f"  Python:   {sys.version}")
    print(f"  Script:   {SCRIPT_DIR}")
    print()

    step(1, 4, "Checking Python version...")
    check_python()

    step(2, 4, "Checking required files...")
    check_files()

    step(3, 4, "Checking Python packages...")
    missing = check_packages()

    step(4, 4, "Checking bioinformatics tools...")
    found, not_found = check_biotools()

    # Summary
    print(f"\n{'='*60}")
    print(f"{C.BOLD}  Environment Summary{C.RESET}")
    print(f"{'='*60}")
    total_tools = len(BIOTOOLS)
    n_found = len(found)
    if not missing and n_found == total_tools:
        ok("All dependencies satisfied — full functionality available")
    elif not missing:
        ok(f"Core dependencies OK — {n_found}/{total_tools} optional tools installed")
        if not_found:
            info("Missing optional tools can be installed via conda:")
            unique_pkgs = set(m["conda"] for m in not_found.values())
            info(f"  conda install -c bioconda -c conda-forge {' '.join(sorted(unique_pkgs))}")
    else:
        warn(f"Missing {len(missing)} Python packages — run: pip install -r requirements.txt")

    print()


# ── Main ─────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description=f"pLIN v{APP_VERSION} — Setup & Launch",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 setup_pLIN.py              Install dependencies + launch pLIN
  python3 setup_pLIN.py --install    Install only (no launch)
  python3 setup_pLIN.py --launch     Launch only (skip install)
  python3 setup_pLIN.py --check      Check environment only
  python3 setup_pLIN.py --docker     Build and run via Docker
        """,
    )
    parser.add_argument("--install", action="store_true", help="Install dependencies only")
    parser.add_argument("--launch", action="store_true", help="Launch pLIN only (skip install)")
    parser.add_argument("--check", action="store_true", help="Check environment only")
    parser.add_argument("--docker", action="store_true", help="Build and run via Docker")
    parser.add_argument("--no-biotools", action="store_true", help="Skip bioinformatics tool installation")
    args = parser.parse_args()

    # Docker mode
    if args.docker:
        docker_build_and_run()
        return

    # Check-only mode
    if args.check:
        environment_report()
        return

    banner()
    plat = detect_platform()
    print(f"  Platform: {plat['os_name']}")
    print(f"  Python:   {sys.version}")
    print()

    if not args.launch:
        # Step 1: Check Python
        step(1, 4, "Checking environment...")
        if not check_python():
            fail(f"Python {MIN_PYTHON[0]}.{MIN_PYTHON[1]}+ is required")
            sys.exit(1)

        if not check_files():
            fail("Required data files are missing. Ensure you have the complete pLIN distribution.")
            sys.exit(1)

        # Step 2: Install Python packages
        install_python_packages()

        # Step 3: Install bioinformatics tools
        if not args.no_biotools:
            install_biotools_conda()
        else:
            step(3, 4, "Skipping bioinformatics tools (--no-biotools)")
            info("Core pLIN features will work; AMR detection requires AMRFinderPlus")

    if not args.install:
        # Step 4: Launch
        launch()
    else:
        print(f"\n{C.BOLD}Installation complete.{C.RESET}")
        print(f"  Launch pLIN with: python3 {__file__} --launch")
        print(f"  Or directly:      streamlit run plin_app.py")
        print()


if __name__ == "__main__":
    main()
