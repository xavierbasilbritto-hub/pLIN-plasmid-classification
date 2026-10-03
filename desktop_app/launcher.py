#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Standalone desktop launcher for the pLIN Streamlit app.

PyInstaller bundles this script (not plin_app.py directly), because a
Streamlit app is normally started via the `streamlit run` CLI rather than
plain `python script.py`. This launcher invokes Streamlit's own internal
CLI programmatically, then opens the user's default browser once the local
server is ready: giving a standalone double-clickable app with the same
"opens in your browser" experience as running `streamlit run plin_app.py`
by hand, but with no separate Python/pip install required by the user.

External bioinformatics tools (AMRFinderPlus, MOB-suite, BLAST+, FastANI,
MinCED, Prodigal, minimap2) are NOT bundled: they are substantial,
platform-specific binaries normally distributed via conda/bioconda, and
plin_app.py already detects their presence at runtime and gracefully
disables the modules that need them if they're not found on the system.
The core pLIN pipeline (4-mer vectorisation, clustering, lineage coding,
replicon classification, outbreak detection, GUI) needs no external tool
and works fully standalone.
"""

import os
import socket
import sys
import threading
import time
import webbrowser

# Streamlit auto-detects "development mode" by checking whether it's being
# run from inside a git checkout of the streamlit package itself; inside a
# PyInstaller bundle this check can misfire, and development mode refuses
# to honour --server.port. Force it off explicitly before Streamlit's own
# config module is imported anywhere (including transitively).
os.environ["STREAMLIT_GLOBAL_DEVELOPMENT_MODE"] = "false"


def resource_path(relative_path):
    """Resolve a path that works both in a PyInstaller onedir/onefile
    bundle (sys._MEIPASS) and when running this script directly."""
    base_path = getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)


def selftest(base_dir):
    """`pLIN --selftest`: check, inside the packaged app, that pLIN v4.1 typing can run (gene prediction
    with pyrodigal, k-mer sketching, the v4.1 modules and the bundled MMseqs2). Exit code 0 when all pass."""
    import subprocess
    sys.path.insert(0, base_dir)
    import pyrodigal
    from Bio import SeqIO
    from plin_kmers import adaptive_sketch
    from plin_v41_typer import PlinV41Release, find_mmseqs  # noqa: F401  (import check)
    fasta = os.path.join(base_dir, "sample_data", "swiss_vim1_outbreak", "NARACHVIM12_plasmids.fasta")
    seq = str(next(SeqIO.parse(fasta, "fasta")).seq)
    genes = pyrodigal.GeneFinder(meta=True).find_genes(seq.encode())
    print(f"pyrodigal: {len(genes)} genes predicted")
    sk, scale = adaptive_sketch(seq)
    print(f"k-mer sketch: {len(sk)} hashes (scale {scale})")
    mm = find_mmseqs()
    print(f"MMseqs2: {mm}")
    if not mm or not len(genes):
        return 1
    r = subprocess.run([mm, "version"], capture_output=True, text=True)
    print(f"MMseqs2 version: {r.stdout.strip()}")
    if r.returncode != 0:
        return 1
    print("SELFTEST PASSED")
    return 0


def find_free_port(preferred=8501):
    """Use the preferred port if free, otherwise let the OS pick one."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        try:
            s.bind(("127.0.0.1", preferred))
            return preferred
        except OSError:
            s.bind(("127.0.0.1", 0))
            return s.getsockname()[1]


def open_browser_when_ready(url, timeout_s=30):
    """Poll the local server until it accepts connections, then open the
    browser. Avoids racing Streamlit's own startup time."""
    deadline = time.time() + timeout_s
    host, port = url.split("://")[1].split(":")
    port = int(port)
    while time.time() < deadline:
        try:
            with socket.create_connection((host, port), timeout=0.5):
                break
        except OSError:
            time.sleep(0.3)
    webbrowser.open(url)


def main():
    # Bundled app lives alongside this launcher inside the PyInstaller
    # bundle root; running un-bundled (e.g. `python launcher.py` during
    # development) falls back to the real repo layout one level up.
    if hasattr(sys, "_MEIPASS"):
        app_path = resource_path("plin_app.py")
        base_dir = resource_path(".")
    else:
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        app_path = os.path.join(base_dir, "plin_app.py")

    # plin_app.py resolves its own data/output paths relative to its own
    # file location and/or the current working directory in several
    # places: run from base_dir so both resolve the same way as a normal
    # `streamlit run plin_app.py` invocation from the repo root would.
    os.chdir(base_dir)

    if "--selftest" in sys.argv:
        sys.exit(selftest(base_dir))

    port = find_free_port(8501)
    url = f"http://localhost:{port}"

    print("=" * 60)
    print("pLIN: plasmid Lineage Identification Number")
    print("=" * 60)
    print(f"Starting local server at {url}")
    print("This window must stay open while pLIN is running.")
    print("Close this window (or press Ctrl+C) to quit pLIN.")
    print("=" * 60)

    threading.Thread(target=open_browser_when_ready, args=(url,), daemon=True).start()

    # Programmatic equivalent of:
    #   streamlit run plin_app.py --server.port <port> --server.headless true
    sys.argv = [
        "streamlit",
        "run",
        app_path,
        "--server.port", str(port),
        "--server.headless", "true",
        "--browser.gatherUsageStats", "false",
        "--server.fileWatcherType", "none",
    ]
    from streamlit.web import cli as stcli
    sys.exit(stcli.main())


if __name__ == "__main__":
    main()
