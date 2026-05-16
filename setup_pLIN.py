#!/usr/bin/env python3
"""
pLIN Tool — Cross-Platform Setup & Launch Script (Single File)
==============================================================
One file to install dependencies and launch pLIN on macOS, Linux, or Windows.

Usage:
    python3 setup_pLIN.py              # Install + launch
    python3 setup_pLIN.py --install    # Install only (no launch)
    python3 setup_pLIN.py --launch     # Launch only (skip install)
    python3 setup_pLIN.py --check      # Check environment only
    python3 setup_pLIN.py --docker     # Build and run Docker container
    python3 setup_pLIN.py --uninstall  # Remove conda env + venv

Works on: macOS (Intel/Apple Silicon), Linux (Ubuntu/Debian/Fedora/Arch/RHEL/openSUSE), Windows 10/11

Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
"""

import os
import sys
import platform
import subprocess
import shutil
import argparse
import webbrowser
import time
import threading
from pathlib import Path

# ── Constants ────────────────────────────────────────────────────────────
APP_NAME = "pLIN"
APP_VERSION = "3.0.0"
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
    "mob_typer":    {"conda": None,                  "desc": "Plasmid mobility typing (pip)"},
    "mlst":         {"conda": "mlst",                "desc": "Multi-locus sequence typing"},
}

# Conda packages to install (deduplicated)
CONDA_PACKAGES = [
    "ncbi-amrfinderplus", "mash", "fastani", "minimap2",
    "minced", "blast", "prodigal", "mlst",
]

# Output directories to create
OUTPUT_DIRS = [
    OUTPUT_DIR / "amrfinder",
    OUTPUT_DIR / "crispr_evaluation",
    OUTPUT_DIR / "figures",
    OUTPUT_DIR / "manuscripts" / "docx",
    OUTPUT_DIR / "mge_detection",
]

# IS element reference accessions for MGE detection database
# 25 representative IS families from ISfinder/NCBI
IS_REFERENCE_ACCESSIONS = {
    # Gram-negative IS families (18)
    "IS26":    "X00011.1",
    "ISEcp1":  "AJ242809.1",
    "IS1":     "J01730.1",
    "IS903":   "M17148.1",
    "IS6100":  "M95400.1",
    "ISKpn26": "KF914891.1",
    "IS5":     "X02311.1",
    "IS3":     "X02180.1",
    "IS4321":  "AJ245418.1",
    "IS15":    "X01840.1",
    "IS10":    "J01830.1",
    "IS2":     "J01733.1",
    "IS4":     "V00029.1",
    "IS30":    "X00792.1",
    "IS66":    "X53365.1",
    "IS110":   "M21395.1",
    "ISPa":    "AF261825.1",
    "ISAba":   "AY758396.1",
    # Gram-positive IS families (7)
    "IS256":   "M18086.1",
    "IS257":   "U40412.1",
    "IS16":    "AF053365.1",
    "ISEnfa":  "AF162694.1",
    "IS1216":  "L40841.1",
    "IS1251":  "X83579.1",
    "Tn916":   "U09422.1",
}

# Miniconda download URLs
MINICONDA_URLS = {
    ("darwin", "arm64"):    "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh",
    ("darwin", "x86_64"):   "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh",
    ("linux", "x86_64"):    "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh",
    ("linux", "aarch64"):   "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh",
    ("linux", "arm64"):     "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh",
    ("windows", "amd64"):   "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe",
    ("windows", "x86_64"):  "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe",
}


# ── Colours for terminal output ──────────────────────────────────────────
class C:
    """ANSI colour codes (disabled on Windows CMD without VT support)."""
    _win_vt = False
    if os.name == "nt":
        try:
            import ctypes
            k = ctypes.windll.kernel32
            k.SetConsoleMode(k.GetStdHandle(-11), 7)
            _win_vt = True
        except Exception:
            pass
    _enabled = sys.stdout.isatty() and (os.name != "nt" or _win_vt)
    BOLD  = "\033[1m"   if _enabled else ""
    GREEN = "\033[92m"  if _enabled else ""
    YELLOW= "\033[93m"  if _enabled else ""
    RED   = "\033[91m"  if _enabled else ""
    CYAN  = "\033[96m"  if _enabled else ""
    RESET = "\033[0m"   if _enabled else ""


def banner():
    print(f"""
{C.CYAN}{'='*62}
  {C.BOLD}pLIN: Plasmid Lineage Identification Number System{C.RESET}{C.CYAN}
  Version {APP_VERSION} — Cross-Platform Setup & Launch
  Hierarchical Plasmid Classification + AMR Surveillance
{'='*62}{C.RESET}
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
    """Detect OS, architecture, and distribution details."""
    system = platform.system().lower()
    machine = platform.machine().lower()
    is_arm = machine in ("arm64", "aarch64")

    result = {
        "system": system,
        "machine": machine,
        "is_arm": is_arm,
        "os_name": "",
        "distro": "",
        "distro_id": "",
        "pkg_mgr": "",
    }

    if system == "darwin":
        result["os_name"] = "macOS (Apple Silicon)" if is_arm else "macOS (Intel)"
        try:
            _, ver = run_cmd(["sw_vers", "-productVersion"])
            result["os_name"] += f" {ver}"
        except Exception:
            pass

    elif system == "linux":
        result["os_name"] = f"Linux ({machine})"
        # Detect distribution
        if Path("/etc/os-release").exists():
            os_release = {}
            for line in Path("/etc/os-release").read_text().splitlines():
                if "=" in line:
                    k, v = line.split("=", 1)
                    os_release[k] = v.strip('"')
            result["distro"] = f"{os_release.get('NAME', 'Linux')} {os_release.get('VERSION_ID', '')}"
            result["distro_id"] = os_release.get("ID", "")
            result["os_name"] = f"{result['distro']} ({machine})"

        # Detect package manager
        for mgr in ["apt-get", "dnf", "yum", "pacman", "zypper"]:
            if shutil.which(mgr):
                result["pkg_mgr"] = mgr.replace("-get", "")
                break

    elif system == "windows":
        result["os_name"] = f"Windows ({machine})"

    return result


# ── Command helpers ──────────────────────────────────────────────────────
def run_cmd(cmd, check=True, capture=True, timeout=600, env=None):
    """Run a shell command and return (success, stdout)."""
    try:
        r = subprocess.run(
            cmd, shell=isinstance(cmd, str), check=check,
            capture_output=capture, text=True, timeout=timeout, env=env,
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

    home = Path.home()
    candidates = [
        home / "miniconda3" / "bin" / "conda",
        home / "miniforge3" / "bin" / "conda",
        home / "anaconda3" / "bin" / "conda",
        home / "mambaforge" / "bin" / "conda",
        home / "opt" / "miniconda3" / "bin" / "conda",
        Path("/opt/conda/bin/conda"),
    ]

    if platform.system().lower() == "darwin":
        candidates += [
            Path("/opt/homebrew/Caskroom/miniconda/base/bin/conda"),
            Path("/usr/local/Caskroom/miniconda/base/bin/conda"),
        ]

    if platform.system() == "Windows":
        candidates += [
            home / "miniconda3" / "Scripts" / "conda.exe",
            home / "Miniconda3" / "Scripts" / "conda.exe",
            home / "anaconda3" / "Scripts" / "conda.exe",
            home / "Anaconda3" / "Scripts" / "conda.exe",
            home / "miniforge3" / "Scripts" / "conda.exe",
            home / "mambaforge" / "Scripts" / "conda.exe",
            Path(os.environ.get("LOCALAPPDATA", "")) / "miniconda3" / "Scripts" / "conda.exe",
            Path(os.environ.get("PROGRAMDATA", "")) / "miniconda3" / "Scripts" / "conda.exe",
        ]

    for c in candidates:
        if c.exists():
            return str(c)
    return None


def find_python():
    """Find a suitable Python 3.10+ executable."""
    for cmd in ["python3.12", "python3.11", "python3.10", "python3", "python"]:
        path = which(cmd)
        if path:
            success, ver = run_cmd([path, "--version"])
            if success:
                parts = ver.split()[-1].split(".")
                try:
                    if (int(parts[0]), int(parts[1])) >= MIN_PYTHON:
                        return path, ver.split()[-1]
                except (ValueError, IndexError):
                    pass
    return None, None


def find_conda_env_python():
    """Find the Python binary inside the conda environment."""
    home = Path.home()
    if platform.system() == "Windows":
        candidates = [
            home / "miniconda3" / "envs" / ENV_NAME / "python.exe",
            home / "Miniconda3" / "envs" / ENV_NAME / "python.exe",
            home / "anaconda3" / "envs" / ENV_NAME / "python.exe",
        ]
    else:
        candidates = [
            home / "miniconda3" / "envs" / ENV_NAME / "bin" / "python3",
            home / "miniforge3" / "envs" / ENV_NAME / "bin" / "python3",
            home / "anaconda3" / "envs" / ENV_NAME / "bin" / "python3",
            home / "mambaforge" / "envs" / ENV_NAME / "bin" / "python3",
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


def check_packages(python_bin=None):
    """Verify Python packages are importable."""
    if python_bin is None:
        python_bin = sys.executable
    missing = []
    for pkg in CORE_PACKAGES:
        import_name = pkg.replace("-", "_")
        if import_name == "scikit_learn":
            import_name = "sklearn"
        if import_name == "python_pptx":
            import_name = "pptx"
        if import_name == "biopython":
            import_name = "Bio"
        success, _ = run_cmd(
            [python_bin, "-c", f"import {import_name}"],
            check=False, capture=True, timeout=10,
        )
        if success:
            ok(f"Python package: {pkg}")
        else:
            warn(f"Missing Python package: {pkg}")
            missing.append(pkg)
    return missing


def check_biotools():
    """Check which bioinformatics tools are available."""
    found = {}
    missing = {}
    for tool, meta in BIOTOOLS.items():
        path = which(tool)
        if not path:
            # Also check inside conda env
            path = _find_tool_in_conda_env(tool)
        if path:
            found[tool] = path
            ok(f"{tool} — {meta['desc']}")
        else:
            missing[tool] = meta
            warn(f"{tool} — {meta['desc']} (not found, optional)")
    return found, missing


def _find_tool_in_conda_env(tool):
    """Check common conda env locations for a tool."""
    home = Path.home()
    if platform.system() == "Windows":
        suffixes = [
            Path("envs") / ENV_NAME / "Scripts" / f"{tool}.exe",
            Path("envs") / ENV_NAME / "bin" / tool,
        ]
        bases = [
            home / "miniconda3", home / "Miniconda3",
            home / "anaconda3", home / "Anaconda3",
        ]
    else:
        suffixes = [Path("envs") / ENV_NAME / "bin" / tool]
        bases = [
            home / "miniconda3", home / "miniforge3",
            home / "anaconda3", home / "mambaforge",
        ]
    for base in bases:
        for suffix in suffixes:
            p = base / suffix
            if p.exists():
                return str(p)
    return None


# ── System dependencies (Linux) ─────────────────────────────────────────
def install_system_deps(plat):
    """Install system-level dependencies on Linux."""
    if plat["system"] != "linux":
        return

    info("Checking system dependencies...")
    pkg_mgr = plat["pkg_mgr"]

    # Check for curl/wget
    if not which("curl") and not which("wget"):
        warn("Neither curl nor wget found. Installing...")
        cmds = {
            "apt": "sudo apt-get update -qq && sudo apt-get install -y -qq curl wget build-essential default-jre-headless",
            "dnf": "sudo dnf install -y curl wget gcc gcc-c++ java-11-openjdk-headless",
            "yum": "sudo yum install -y curl wget gcc gcc-c++ java-11-openjdk-headless",
            "pacman": "sudo pacman -Sy --noconfirm curl wget base-devel jre-openjdk-headless",
            "zypper": "sudo zypper install -y curl wget gcc gcc-c++ java-11-openjdk-headless",
        }
        if pkg_mgr in cmds:
            run_cmd(cmds[pkg_mgr], check=False, capture=False)

    # Java is needed for MinCED
    if not which("java"):
        warn("Java not found (needed for MinCED). Installing...")
        cmds = {
            "apt": "sudo apt-get install -y -qq default-jre-headless",
            "dnf": "sudo dnf install -y java-11-openjdk-headless",
            "yum": "sudo yum install -y java-11-openjdk-headless",
            "pacman": "sudo pacman -Sy --noconfirm jre-openjdk-headless",
            "zypper": "sudo zypper install -y java-11-openjdk-headless",
        }
        if pkg_mgr in cmds:
            run_cmd(cmds[pkg_mgr], check=False, capture=False)
        else:
            warn("Install Java manually for MinCED support")

    ok("System dependencies checked")


# ── Installation ─────────────────────────────────────────────────────────
def install_conda_if_missing(plat):
    """Download and install Miniconda if conda is not found."""
    system = plat["system"]
    machine = plat["machine"]

    # Normalise machine name for URL lookup
    url_machine = machine
    if machine == "arm64":
        url_machine = "arm64" if system == "darwin" else "aarch64"
    if machine in ("amd64", "x86_64", "x64"):
        url_machine = "x86_64"

    key = (system, url_machine)
    if system == "windows":
        key = ("windows", "x86_64")

    url = MINICONDA_URLS.get(key)
    if not url:
        fail(f"No Miniconda installer available for {system}/{machine}")
        return None

    info(f"Conda not found. Downloading Miniconda for {system}/{machine}...")

    if system == "windows":
        installer = Path.home() / "miniconda_installer.exe"
        success, _ = run_cmd(
            f'powershell -Command "Invoke-WebRequest -Uri {url} -OutFile {installer}"',
            check=False, capture=True, timeout=300,
        )
        if success and installer.exists():
            info("Running Miniconda installer (silent)...")
            run_cmd(
                f'start /wait "" "{installer}" /InstallationType=JustMe /RegisterPython=0 '
                f'/AddToPath=0 /S /D={Path.home() / "miniconda3"}',
                check=False, capture=False, timeout=600,
            )
            installer.unlink(missing_ok=True)
            conda_bin = str(Path.home() / "miniconda3" / "Scripts" / "conda.exe")
        else:
            fail("Failed to download Miniconda")
            return None
    else:
        installer = Path("/tmp/miniconda_installer.sh")
        dl_cmd = f'curl -fsSL "{url}" -o "{installer}"'
        if not which("curl"):
            dl_cmd = f'wget -q "{url}" -O "{installer}"'
        success, _ = run_cmd(dl_cmd, check=False, capture=True, timeout=300)
        if success and installer.exists():
            run_cmd(f'bash "{installer}" -b -p {Path.home() / "miniconda3"}',
                    check=False, capture=False, timeout=300)
            installer.unlink(missing_ok=True)
            conda_bin = str(Path.home() / "miniconda3" / "bin" / "conda")
        else:
            fail("Failed to download Miniconda")
            return None

    if Path(conda_bin).exists():
        ok(f"Miniconda installed: {conda_bin}")
        return conda_bin
    else:
        fail("Miniconda installation failed")
        return None


def setup_conda_env(conda):
    """Create the pLIN_tools conda environment if it doesn't exist."""
    # Check if env exists
    success, envs = run_cmd([conda, "env", "list"], check=False)
    if success and ENV_NAME in envs:
        info(f"Conda environment '{ENV_NAME}' already exists")
    else:
        info(f"Creating conda environment '{ENV_NAME}' with Python 3.11...")
        success, out = run_cmd(
            [conda, "create", "-n", ENV_NAME, "python=3.11", "-y", "--quiet"],
            check=False, capture=True, timeout=300,
        )
        if success:
            ok("Conda environment created")
        else:
            fail(f"Failed to create conda env: {out[:200]}")
            return False

    # Configure channels
    for ch in ["defaults", "bioconda", "conda-forge"]:
        run_cmd([conda, "config", "--add", "channels", ch], check=False)
    run_cmd([conda, "config", "--set", "channel_priority", "strict"], check=False)

    return True


def install_python_packages(python_bin=None):
    """Install Python dependencies via pip."""
    if python_bin is None:
        python_bin = sys.executable

    info("Upgrading pip...")
    run_cmd([python_bin, "-m", "pip", "install", "--upgrade", "pip", "--quiet"],
            check=False, capture=True)

    info("Installing Python packages from requirements.txt...")
    success, output = run_cmd(
        [python_bin, "-m", "pip", "install", "--quiet", "-r", str(REQUIREMENTS)],
        check=False, capture=True, timeout=300,
    )
    if success:
        ok("All Python packages installed")
        return True
    else:
        fail(f"pip install failed: {output[:200]}")
        info("Try manually: pip install -r requirements.txt")
        return False


def install_biotools_conda(conda):
    """Install bioinformatics tools via conda into the pLIN_tools env."""
    installed = 0
    total = len(CONDA_PACKAGES) + 1  # +1 for MOBsuite

    for pkg in CONDA_PACKAGES:
        info(f"Installing {pkg}...")
        success, _ = run_cmd(
            [conda, "install", "-n", ENV_NAME, "-y", "-c", "bioconda",
             "-c", "conda-forge", pkg, "--quiet"],
            check=False, capture=True, timeout=300,
        )
        if success:
            ok(f"Installed {pkg}")
            installed += 1
        else:
            warn(f"Could not install {pkg} (optional)")

    # MOBsuite via pip (not available in conda)
    info("Installing MOBsuite via pip...")
    env_python = find_conda_env_python()
    if env_python:
        success, _ = run_cmd(
            [env_python, "-m", "pip", "install", "mob_suite", "--quiet"],
            check=False, capture=True, timeout=300,
        )
        if success:
            ok("Installed MOBsuite")
            installed += 1
        else:
            warn("Could not install MOBsuite (optional)")
    else:
        warn("Could not find conda env Python for MOBsuite install")

    # Update AMRFinderPlus database
    amrfinder = which("amrfinder") or _find_tool_in_conda_env("amrfinder")
    if amrfinder:
        info("Updating AMRFinderPlus database...")
        run_cmd([amrfinder, "--update"], check=False, capture=True, timeout=120)

    ok(f"Installed {installed}/{total} bioinformatics tools")
    return installed


def setup_directories():
    """Create output directory structure."""
    for d in OUTPUT_DIRS:
        d.mkdir(parents=True, exist_ok=True)
    ok("Output directories ready")


def setup_is_reference_database(python_bin=None):
    """Download IS element reference sequences and build BLAST database.

    Downloads curated IS element sequences from NCBI nucleotide database
    using BioPython's Entrez, writes them as a combined FASTA file, and
    builds a BLAST nucleotide database for IS element detection.

    This is optional — if NCBI is unreachable or BioPython is unavailable,
    the function logs a warning and returns False.
    """
    mge_dir = OUTPUT_DIR / "mge_detection"
    mge_dir.mkdir(parents=True, exist_ok=True)

    is_fasta = mge_dir / "is_reference_sequences.fasta"
    is_db_flag = mge_dir / "is_reference_sequences.fasta.ndb"

    # Skip if database already exists
    if is_fasta.exists() and is_fasta.stat().st_size > 1000 and is_db_flag.exists():
        ok("IS reference database already exists")
        return True

    if python_bin is None:
        python_bin = sys.executable

    info("Downloading IS element reference sequences from NCBI...")
    info("(Requires internet — skipped if unavailable)")

    # Inline download script to keep setup_pLIN.py self-contained
    download_script = (
        "import sys\n"
        "try:\n"
        "    from Bio import Entrez, SeqIO\n"
        "except ImportError:\n"
        '    print("SKIP:BioPython not available")\n'
        "    sys.exit(0)\n"
        'Entrez.email = "plin_tool@example.com"\n'
        f"accessions = {repr(IS_REFERENCE_ACCESSIONS)}\n"
        f'out_path = "{str(is_fasta)}"\n'
        "records = []\n"
        "for is_name, acc in accessions.items():\n"
        "    try:\n"
        '        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")\n'
        '        rec = SeqIO.read(handle, "fasta")\n'
        "        handle.close()\n"
        "        rec.id = is_name\n"
        "        rec.description = is_name\n"
        "        records.append(rec)\n"
        "    except Exception as e:\n"
        '        print(f"WARN:{is_name}: {e}")\n'
        "if records:\n"
        '    with open(out_path, "w") as f:\n'
        '        SeqIO.write(records, f, "fasta")\n'
        '    print(f"OK:{len(records)}/{len(accessions)} IS references downloaded")\n'
        "else:\n"
        '    print("FAIL:No sequences downloaded")\n'
        "    sys.exit(1)\n"
    )

    success, output = run_cmd(
        [python_bin, "-c", download_script],
        check=False, capture=True, timeout=120,
    )

    if not success or "FAIL:" in (output or ""):
        warn("Could not download IS reference sequences (NCBI may be unreachable)")
        warn("IS element detection can be set up later with: python3 detect_IS_elements.py")
        return False

    if "SKIP:" in (output or ""):
        warn("BioPython not available — skipping IS reference database setup")
        return False

    ok(output.split("OK:")[-1] if "OK:" in output else "IS references downloaded")

    # Build BLAST database
    makeblastdb_bin = which("makeblastdb") or _find_tool_in_conda_env("makeblastdb")
    if makeblastdb_bin and is_fasta.exists():
        info("Building BLAST database for IS element detection...")
        success, output = run_cmd(
            [makeblastdb_bin, "-in", str(is_fasta), "-dbtype", "nucl",
             "-parse_seqids", "-out", str(is_fasta)],
            check=False, capture=True, timeout=60,
        )
        if success:
            ok("IS element BLAST database built")
            return True
        else:
            warn(f"makeblastdb failed: {(output or '')[:200]}")
            return False
    else:
        warn("makeblastdb not found — BLAST database not built")
        warn("Install BLAST and run: makeblastdb -in output/mge_detection/is_reference_sequences.fasta -dbtype nucl")
        return False


def verify_blast_for_contig_classification():
    """Verify BLAST tools are available for IS detection and contig analysis.

    The contig classifier itself uses 4-mer KNN (no BLAST needed), but
    IS element detection and CRISPR host prediction both require BLAST.
    """
    blastn_bin = which("blastn") or _find_tool_in_conda_env("blastn")
    makeblastdb_bin = which("makeblastdb") or _find_tool_in_conda_env("makeblastdb")

    if blastn_bin and makeblastdb_bin:
        success, version = run_cmd([blastn_bin, "-version"], check=False, timeout=10)
        if success:
            ver_line = version.split("\n")[0] if version else "unknown"
            ok(f"BLAST available: {ver_line}")
            return True

    warn("BLAST tools not found — IS element detection will be unavailable")
    info("Core pLIN classification and contig filtering work without BLAST")
    return False


# ── Docker ───────────────────────────────────────────────────────────────
def docker_build_and_run():
    """Build and run pLIN via Docker."""
    banner()
    step(1, 3, "Checking Docker...")

    if not which("docker"):
        fail("Docker not found")
        system = platform.system().lower()
        if system == "linux":
            info("Install Docker:")
            info("  Ubuntu/Debian: sudo apt install docker.io docker-compose")
            info("  Fedora/RHEL:   sudo dnf install docker docker-compose")
            info("  Arch:          sudo pacman -S docker docker-compose")
            info("  Or:            https://docs.docker.com/engine/install/")
        else:
            info("Install Docker Desktop from: https://www.docker.com/products/docker-desktop/")
        sys.exit(1)

    ok("Docker found")

    step(2, 3, "Building pLIN Docker image (this may take 5-10 minutes on first run)...")
    success, output = run_cmd(
        f"docker build -t plin:latest {SCRIPT_DIR}",
        check=False, capture=False, timeout=600,
    )
    if not success:
        fail("Docker build failed")
        sys.exit(1)
    ok("Docker image built: plin:latest")

    step(3, 3, "Starting pLIN container...")
    run_cmd("docker stop plin-app 2>/dev/null", check=False)
    run_cmd("docker rm plin-app 2>/dev/null", check=False)

    vol_data = f"{SCRIPT_DIR / 'data'}:/app/data"
    vol_output = f"{SCRIPT_DIR / 'output'}:/app/output"
    success, output = run_cmd(
        f'docker run -d --name plin-app -p {STREAMLIT_PORT}:{STREAMLIT_PORT} '
        f'-v "{vol_data}" -v "{vol_output}" plin:latest',
        check=False, capture=True,
    )
    if success:
        ok("pLIN running in Docker container")
        print(f"\n{C.BOLD}  Open in browser: http://localhost:{STREAMLIT_PORT}{C.RESET}\n")
        info("Stop with:  docker stop plin-app")
        info("Logs:       docker logs -f plin-app")
        _open_browser()
    else:
        fail(f"Failed to start container: {output[:200]}")


# ── Uninstall ────────────────────────────────────────────────────────────
def uninstall():
    """Remove conda environment and virtual environment."""
    banner()
    print(f"  {C.YELLOW}This will remove the pLIN conda environment and virtual environment.{C.RESET}")
    print(f"  {C.YELLOW}Your data files and sequences will NOT be deleted.{C.RESET}")
    print()

    try:
        confirm = input("  Continue? [y/N]: ").strip().lower()
    except (EOFError, KeyboardInterrupt):
        print("\n  Cancelled.")
        return

    if confirm != "y":
        print("  Cancelled.")
        return

    # Remove conda env
    conda = find_conda()
    if conda:
        success, envs = run_cmd([conda, "env", "list"], check=False)
        if success and ENV_NAME in envs:
            info(f"Removing conda environment '{ENV_NAME}'...")
            run_cmd([conda, "env", "remove", "-n", ENV_NAME, "-y"], check=False)
            ok("Conda environment removed")

    # Remove venv
    venv_dir = SCRIPT_DIR / ".venv"
    if venv_dir.exists():
        info("Removing virtual environment...")
        shutil.rmtree(venv_dir)
        ok("Virtual environment removed")

    ok("Uninstall complete. Data files are preserved.")


# ── Launch ───────────────────────────────────────────────────────────────
def _open_browser():
    """Open the browser after a short delay."""
    def _open():
        time.sleep(3)
        try:
            webbrowser.open(f"http://localhost:{STREAMLIT_PORT}")
        except Exception:
            pass
    t = threading.Thread(target=_open, daemon=True)
    t.start()


def _find_streamlit():
    """Find the streamlit executable."""
    if which("streamlit"):
        return which("streamlit")

    # Check conda env
    home = Path.home()
    if platform.system() == "Windows":
        candidates = [
            home / "miniconda3" / "envs" / ENV_NAME / "Scripts" / "streamlit.exe",
            home / "Miniconda3" / "envs" / ENV_NAME / "Scripts" / "streamlit.exe",
            home / "anaconda3" / "envs" / ENV_NAME / "Scripts" / "streamlit.exe",
        ]
    else:
        candidates = [
            home / "miniconda3" / "envs" / ENV_NAME / "bin" / "streamlit",
            home / "miniforge3" / "envs" / ENV_NAME / "bin" / "streamlit",
            home / "anaconda3" / "envs" / ENV_NAME / "bin" / "streamlit",
            home / "mambaforge" / "envs" / ENV_NAME / "bin" / "streamlit",
        ]

    # Check venv
    if platform.system() == "Windows":
        candidates.append(SCRIPT_DIR / ".venv" / "Scripts" / "streamlit.exe")
    else:
        candidates.append(SCRIPT_DIR / ".venv" / "bin" / "streamlit")

    for c in candidates:
        if c.exists():
            return str(c)
    return None


def launch():
    """Launch the Streamlit app."""
    if not APP_FILE.exists():
        fail(f"Application file not found: {APP_FILE}")
        sys.exit(1)

    streamlit = _find_streamlit()
    if not streamlit:
        fail("Streamlit not found. Run: python3 setup_pLIN.py --install")
        sys.exit(1)

    print(f"""
{C.BOLD}  ╔══════════════════════════════════════════════════╗
  ║  pLIN is starting at: http://localhost:{STREAMLIT_PORT}    ║
  ║  Press Ctrl+C to stop the server.               ║
  ╚══════════════════════════════════════════════════╝{C.RESET}
""")

    _open_browser()

    try:
        subprocess.run(
            [
                streamlit, "run", str(APP_FILE),
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
    if plat["distro"]:
        print(f"  Distro:   {plat['distro']}")
    if plat["pkg_mgr"]:
        print(f"  Pkg mgr:  {plat['pkg_mgr']}")
    print(f"  Python:   {sys.version}")
    print(f"  Script:   {SCRIPT_DIR}")
    print()

    step(1, 6, "Checking Python version...")
    check_python()

    step(2, 6, "Checking required files...")
    check_files()

    step(3, 6, "Checking Python packages...")
    # Try conda env Python first, fall back to system Python
    env_python = find_conda_env_python()
    if env_python:
        info(f"Using conda env Python: {env_python}")
        missing = check_packages(env_python)
    else:
        missing = check_packages()

    step(4, 6, "Checking bioinformatics tools...")
    found, not_found = check_biotools()

    step(5, 6, "Checking conda...")
    conda = find_conda()
    if conda:
        ok(f"Conda: {conda}")
        success, envs = run_cmd([conda, "env", "list"], check=False)
        if success and ENV_NAME in envs:
            ok(f"Conda environment '{ENV_NAME}' exists")
        else:
            warn(f"Conda environment '{ENV_NAME}' not found")
    else:
        warn("Conda not found")

    step(6, 6, "Checking IS reference database...")
    is_fasta = OUTPUT_DIR / "mge_detection" / "is_reference_sequences.fasta"
    is_db_flag = OUTPUT_DIR / "mge_detection" / "is_reference_sequences.fasta.ndb"
    if is_fasta.exists() and is_db_flag.exists():
        ok("IS element reference database found")
    elif is_fasta.exists():
        warn("IS reference FASTA exists but BLAST database not built")
        info("Run: makeblastdb -in output/mge_detection/is_reference_sequences.fasta -dbtype nucl")
    else:
        warn("IS element reference database not set up")
        info("Run: python3 setup_pLIN.py --install")
    verify_blast_for_contig_classification()

    # Summary
    print(f"\n{'='*62}")
    print(f"{C.BOLD}  Environment Summary{C.RESET}")
    print(f"{'='*62}")
    total_tools = len(BIOTOOLS)
    n_found = len(found)
    if not missing and n_found == total_tools:
        ok("All dependencies satisfied — full functionality available")
    elif not missing:
        ok(f"Core dependencies OK — {n_found}/{total_tools} optional tools installed")
        if not_found:
            unique_pkgs = set()
            for m in not_found.values():
                if m["conda"]:
                    unique_pkgs.add(m["conda"])
            if unique_pkgs:
                info("Missing optional tools can be installed via conda:")
                info(f"  conda install -n {ENV_NAME} -c bioconda -c conda-forge {' '.join(sorted(unique_pkgs))}")
            mob_missing = "mob_typer" in not_found
            if mob_missing:
                info("Install MOBsuite: pip install mob_suite")
    else:
        warn(f"Missing {len(missing)} Python packages — run: python3 setup_pLIN.py --install")

    print()


# ── Main ─────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description=f"pLIN v{APP_VERSION} — Cross-Platform Setup & Launch (macOS / Linux / Windows)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Examples:
  python3 setup_pLIN.py              Install dependencies + launch pLIN
  python3 setup_pLIN.py --install    Install only (no launch)
  python3 setup_pLIN.py --launch     Launch only (skip install)
  python3 setup_pLIN.py --check      Check environment only
  python3 setup_pLIN.py --docker     Build and run via Docker
  python3 setup_pLIN.py --uninstall  Remove conda env + venv

Supported platforms:
  macOS     Ventura 13+, Sonoma 14+, Sequoia 15+ (Intel & Apple Silicon)
  Linux     Ubuntu 20.04+, Debian 11+, Fedora 36+, CentOS/RHEL 8+,
            Arch Linux, openSUSE 15+ (x86_64 and ARM64)
  Windows   Windows 10 (1903+), Windows 11 (x64)
            For full bioinformatics support, use WSL2: wsl --install
""",
    )
    parser.add_argument("--install", action="store_true", help="Install dependencies only")
    parser.add_argument("--launch", action="store_true", help="Launch pLIN only (skip install)")
    parser.add_argument("--check", action="store_true", help="Check environment only")
    parser.add_argument("--docker", action="store_true", help="Build and run via Docker")
    parser.add_argument("--uninstall", action="store_true", help="Remove conda env + venv")
    parser.add_argument("--no-biotools", action="store_true", help="Skip bioinformatics tool installation")
    args = parser.parse_args()

    # Docker mode
    if args.docker:
        docker_build_and_run()
        return

    # Uninstall mode
    if args.uninstall:
        uninstall()
        return

    # Check-only mode
    if args.check:
        environment_report()
        return

    # Launch-only mode
    if args.launch:
        banner()
        step(1, 1, "Launching pLIN...")
        launch()
        return

    # ── Full install (or install-only) ───────────────────────────────────
    banner()
    plat = detect_platform()
    total_steps = 7
    print(f"  Platform: {plat['os_name']}")
    print(f"  Python:   {sys.version}")
    print()

    # Step 1: Check environment
    step(1, total_steps, "Checking environment...")
    if not check_python():
        warn(f"Python {MIN_PYTHON[0]}.{MIN_PYTHON[1]}+ required — will install via conda")

    if not check_files():
        fail("Required data files are missing. Ensure you have the complete pLIN distribution.")
        sys.exit(1)

    # Step 2: System dependencies (Linux only)
    step(2, total_steps, "Checking system dependencies...")
    install_system_deps(plat)
    if plat["system"] != "linux":
        ok("No system dependencies needed")

    # Step 3: Conda environment + bioinformatics tools
    step(3, total_steps, "Setting up conda environment + bioinformatics tools...")
    if not args.no_biotools:
        conda = find_conda()
        if not conda:
            conda = install_conda_if_missing(plat)
        if conda:
            ok(f"Conda: {conda}")
            if setup_conda_env(conda):
                install_biotools_conda(conda)
        else:
            warn("Conda not available — skipping bioinformatics tools")
            info("Core pLIN features (classification, pLIN assignment) work without these tools")
            info("Install conda from: https://docs.conda.io/en/latest/miniconda.html")
            if plat["system"] == "windows":
                info("For full bioinformatics support on Windows, consider WSL2:")
                info("  1. wsl --install")
                info("  2. Run this script inside WSL")
    else:
        info("Skipping bioinformatics tools (--no-biotools)")
        info("Core pLIN features will work; AMR detection requires AMRFinderPlus")

    # Step 4: Install Python packages
    step(4, total_steps, "Installing Python packages...")
    env_python = find_conda_env_python()
    if env_python:
        info(f"Using conda env Python: {env_python}")
        install_python_packages(env_python)
    else:
        install_python_packages()

    # Step 5: Setup directories
    step(5, total_steps, "Setting up directories...")
    setup_directories()

    # Step 6: IS reference database + BLAST verification
    step(6, total_steps, "Setting up IS reference database + verifying BLAST...")
    env_python = find_conda_env_python() or sys.executable
    if not args.no_biotools:
        setup_is_reference_database(env_python)
        verify_blast_for_contig_classification()
    else:
        info("Skipping IS database setup (--no-biotools)")

    # Step 7: Launch or finish
    if args.install:
        print(f"\n{C.BOLD}{'='*62}{C.RESET}")
        print(f"{C.BOLD}  Installation Complete!{C.RESET}")
        print(f"{C.BOLD}{'='*62}{C.RESET}")
        print(f"  Launch pLIN with: {C.CYAN}python3 {__file__} --launch{C.RESET}")
        print(f"  Or directly:      {C.CYAN}streamlit run plin_app.py{C.RESET}")
        if find_conda():
            print(f"  Or with conda:    {C.CYAN}conda activate {ENV_NAME} && streamlit run plin_app.py{C.RESET}")
        if plat["system"] == "windows":
            print(f"\n  For full bioinformatics support, use WSL2:")
            print(f"    wsl --install")
        print()
    else:
        step(7, total_steps, "Launching pLIN...")
        launch()


if __name__ == "__main__":
    main()
