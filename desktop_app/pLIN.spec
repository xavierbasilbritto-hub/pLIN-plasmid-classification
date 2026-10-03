# -*- mode: python ; coding: utf-8 -*-
# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
#
# PyInstaller spec for the standalone pLIN desktop app.
#
# Streamlit apps need two things PyInstaller does not auto-detect:
#   1. Streamlit's own static frontend assets (the compiled React UI),
#      via collect_data_files("streamlit").
#   2. Package metadata (.dist-info) for streamlit and a few libraries
#      that introspect their own version at runtime, via copy_metadata().
# Both are collected explicitly below; omitting either is the single most
# common cause of a PyInstaller-bundled Streamlit app failing at startup
# or serving a blank page.
#
# Build (run from the PLASMID_TOOL repo root, not from desktop_app/):
#   pyinstaller desktop_app/pLIN.spec --noconfirm

import glob
import os
import re
from PyInstaller.utils.hooks import collect_data_files, copy_metadata, collect_submodules

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(SPEC)), ".."))

block_cipher = None


def _read_plin_app_version():
    """Read PLIN_APP_VERSION from plin_app.py's source text via regex,
    rather than importing the module (which would pull in Streamlit,
    pandas, etc. into the PyInstaller spec's own execution context).
    Single source of truth: bump PLIN_APP_VERSION in plin_app.py, and both
    the app's own provenance stamps (see get_run_provenance()) and this
    bundle's CFBundleShortVersionString stay in sync automatically,
    instead of two hardcoded strings drifting apart (as happened between
    v3.2.0 and v3.2.1, where this string was left stale after a release).
    """
    plin_app_path = os.path.join(REPO_ROOT, "plin_app.py")
    with open(plin_app_path, "r") as f:
        for line in f:
            m = re.match(r'^PLIN_APP_VERSION\s*=\s*"([^"]+)"', line)
            if m:
                return m.group(1)
    raise RuntimeError("Could not find PLIN_APP_VERSION in plin_app.py")


PLIN_APP_VERSION = _read_plin_app_version()

binaries = []
datas = []

# Bundled AMRFinderPlus (macOS/Linux builds only: see build-desktop-app.yml
# for why Windows is excluded: bioconda has no native win-64 build). The CI
# workflow installs AMRFinderPlus into its own conda env before this spec
# runs and exports that env's root as PLIN_AMRFINDER_ENV_DIR; locally, a
# developer can set the same variable (e.g. to their own conda env root,
# such as ~/miniforge3) to build a bundle that includes it. If the variable
# is unset, this app builds exactly as before (AMRFinderPlus detected on
# the user's own system at runtime, per plin_app.py's detect_amrfinder()).
amrfinder_env_dir = os.environ.get("PLIN_AMRFINDER_ENV_DIR")
if amrfinder_env_dir and os.path.isdir(amrfinder_env_dir):
    bundled_bin_dest = "amrfinder_bin"
    # amrfinder itself, its own internal helper binaries (amr_report,
    # fasta_check, etc.: confirmed by running a real end-to-end scan
    # against a bundled build and observing exactly which missing binary
    # it shelled out for next; there is no documented complete list), and
    # the external blastn/blastp/blastx/hmmsearch it also depends on, see
    # plin_app.py's run_amrfinder_on_files(), which passes
    # --blast_bin/--hmmer_bin explicitly at this same bundled path rather
    # than relying on PATH, since a packaged app should not depend on the
    # end user's PATH containing anything.
    for tool in (
        "amrfinder", "amr_report", "amrfinder_index", "amrfinder_update",
        "dna_mutation", "fasta2parts", "fasta_check", "gff_check",
        "blastn", "blastp", "blastx", "tblastn", "hmmsearch",
    ):
        tool_path = os.path.join(amrfinder_env_dir, "bin", tool)
        if os.path.isfile(tool_path):
            binaries.append((tool_path, bundled_bin_dest))
        else:
            print(f"WARNING: bundled-AMRFinderPlus tool not found, skipping: {tool_path}")

    # AMRFinderPlus's gene database (~242MB): a dated subdirectory of
    # share/amrfinderplus/data/. Bundle the lexicographically-latest one,
    # matching plin_app.py's own detect_amrfinder() version-selection logic.
    db_root = os.path.join(amrfinder_env_dir, "share", "amrfinderplus", "data")
    db_versions = sorted(
        d for d in glob.glob(os.path.join(db_root, "*"))
        if os.path.isdir(d) and os.path.basename(d).startswith("20")
    )
    if db_versions:
        latest_db = db_versions[-1]
        datas.append((latest_db, os.path.join("amrfinder_db", os.path.basename(latest_db))))
    else:
        print(f"WARNING: no AMRFinderPlus database version found under {db_root}")
else:
    print("PLIN_AMRFINDER_ENV_DIR not set: building without a bundled AMRFinderPlus "
          "(app will fall back to detecting a system install at runtime, as before).")
# Bundled MMseqs2 (all platforms), used by pLIN v4.1 to place proteins that are new to the
# database. The CI workflow downloads the official static build of MMseqs2 release 18-8cc5c (the
# version that built the database) and exports its unpacked "mmseqs" directory (containing bin/) as
# PLIN_MMSEQS_DIR; plin_v41_typer.find_mmseqs() looks for it at <bundle>/mmseqs/bin/.
mmseqs_dir = os.environ.get("PLIN_MMSEQS_DIR")
if mmseqs_dir and os.path.isdir(os.path.join(mmseqs_dir, "bin")):
    for f in sorted(os.listdir(os.path.join(mmseqs_dir, "bin"))):
        full = os.path.join(mmseqs_dir, "bin", f)
        if os.path.isfile(full):
            binaries.append((full, os.path.join("mmseqs", "bin")))
    if os.path.isfile(os.path.join(mmseqs_dir, "mmseqs.bat")):         # Windows launcher (sets up BusyBox)
        datas.append((os.path.join(mmseqs_dir, "mmseqs.bat"), "mmseqs"))
else:
    print("PLIN_MMSEQS_DIR not set: building without a bundled MMseqs2 (v4.1 typing of plasmids with new "
          "proteins then needs MMseqs2 on the user's system).")

datas += collect_data_files("streamlit")
datas += copy_metadata("streamlit")
datas += copy_metadata("altair")  # streamlit's charting dep also introspects its own metadata
datas += copy_metadata("pandas")
datas += copy_metadata("numpy")

# Bundle the actual pLIN application + the data it ships with, at the
# bundle root so launcher.py's resource_path()/chdir() logic finds them
# exactly where a normal `streamlit run plin_app.py` from the repo root
# would.
datas += [
    (os.path.join(REPO_ROOT, "plin_app.py"), "."),
    # pLIN founder-assignment module imported by plin_app.py at runtime, and
    # the release founder trees query mode places new plasmids against.
    (os.path.join(REPO_ROOT, "plin_founder.py"), "."),
    (os.path.join(REPO_ROOT, "plin_backbone.py"), "."),
    # pLIN v4.1 typing (plin_v41_typer) and the modules it imports. The v4.1 database itself (1.6 GB)
    # is not bundled: the app downloads it on first use into ~/.plin/plin_v41.
    (os.path.join(REPO_ROOT, "plin_kmers.py"), "."),
    (os.path.join(REPO_ROOT, "plin_v4.py"), "."),
    (os.path.join(REPO_ROOT, "plin_v41.py"), "."),
    (os.path.join(REPO_ROOT, "plin_v41_nn.py"), "."),
    (os.path.join(REPO_ROOT, "plin_v41_typer.py"), "."),
    (os.path.join(REPO_ROOT, "validate_alignment_backbone.py"), "."),
    (os.path.join(REPO_ROOT, "data", "plin_founder_tree_training.npz"), "data"),
    (os.path.join(REPO_ROOT, "data", "plin_founder_tree_reference.npz"), "data"),
    (os.path.join(REPO_ROOT, "data", "inc_classifier.npz"), "data"),
    (os.path.join(REPO_ROOT, "data", "inc_centroids.npz"), "data"),
    (os.path.join(REPO_ROOT, "output", "pLIN_assignments.tsv"), "output"),
    # Full expanded reference database, so "query mode"
    # (looking up a newly-uploaded plasmid's nearest neighbour) runs
    # against the full database rather than just the 8,077-plasmid
    # training set: see _load_reference_for_query() in plin_app.py,
    # which already prefers these two files over the smaller training-set
    # fallback when both are present.
    (os.path.join(REPO_ROOT, "output", "pLIN_reference_assignments.tsv"), "output"),
    (os.path.join(REPO_ROOT, "output", "reference_kmer_vectors.npz"), "output"),
    (os.path.join(REPO_ROOT, "assets"), "assets"),
    (os.path.join(REPO_ROOT, ".streamlit", "config.toml"), ".streamlit"),
    # Ready-to-run example dataset (the Swiss VIM-1 outbreak cluster from
    # the manuscript's case study, ~3MB) so a first-time user can try the
    # app immediately without needing their own FASTA files.
    (os.path.join(REPO_ROOT, "sample_data", "swiss_vim1_outbreak"), "sample_data/swiss_vim1_outbreak"),
]

hiddenimports = []
hiddenimports += collect_submodules("streamlit")
hiddenimports += collect_submodules("Bio")
hiddenimports += ["pyrodigal"]   # backbone protein comparison (plin_backbone.py)
hiddenimports += ["scipy.stats", "requests"]   # imported by the v4.1 modules bundled as data files above
hiddenimports += collect_submodules("pyrodigal")   # gene prediction for v4.1 typing
datas += collect_data_files("pyrodigal")
# pyrodigal loads CPU-specific compiled backends from pyrodigal/impl/ (a namespace package without
# __init__.py, so collect_submodules misses them): bundle every compiled module found there.
import pyrodigal as _pyrodigal
_impl = os.path.join(os.path.dirname(_pyrodigal.__file__), "impl")
for _f in glob.glob(os.path.join(_impl, "*.so")) + glob.glob(os.path.join(_impl, "*.pyd")):
    binaries.append((_f, os.path.join("pyrodigal", "impl")))
    hiddenimports.append("pyrodigal.impl." + os.path.basename(_f).split(".")[0])
hiddenimports += [
    "sklearn.utils._typedefs",
    "sklearn.neighbors._partition_nodes",
    "sklearn.utils._heap",
    "sklearn.utils._sorting",
    "sklearn.utils._vector_sentinel",
    "Bio.SeqIO",
    "Bio.Seq",
    "Bio.SeqRecord",
    # matplotlib picks its rendering backend at runtime via
    # importlib.import_module() (e.g. fig.savefig(..., format="pdf")),
    # which PyInstaller's static analysis cannot see. plin_app.py's
    # fig_to_bytes() saves figures as both "png" and "pdf", so both
    # backends must be declared explicitly or the PDF path raises
    # ModuleNotFoundError: No module named 'matplotlib.backends.backend_pdf'
    # the first time a user downloads a PDF or the "download all results"
    # ZIP (which also crashes mid-export, before finishing the AMR figures).
    "matplotlib.backends.backend_agg",
    "matplotlib.backends.backend_pdf",
    "matplotlib.backends.backend_svg",
    "matplotlib.backends.backend_ps",
]
datas += collect_data_files("Bio")

a = Analysis(
    ["launcher.py"],
    pathex=[os.path.dirname(os.path.abspath(SPEC))],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    # torch/transformers back the optional Nucleotide Transformer (LLM)
    # feature, which plin_app.py already imports lazily/on-demand inside
    # the specific function that needs it, not at module load time. They
    # are excluded here because they add ~500MB+ to the bundle for a
    # feature outside this standalone build's scope (see desktop_app
    # build notes); a user who clicks that specific optional feature in
    # the standalone app will see an import error rather than the feature
    # silently working, exactly as intended for this cut-down build.
    excludes=["tkinter", "torch", "torchvision", "torchaudio", "transformers"],
    noarchive=False,
    cipher=block_cipher,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="pLIN",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=True,  # keep a terminal window so the user can see server status / Ctrl+C to quit
    icon=os.path.join(REPO_ROOT, "assets", "pLIN_favicon.png") if os.path.exists(
        os.path.join(REPO_ROOT, "assets", "pLIN_favicon.png")) else None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    name="pLIN",
)

app = BUNDLE(
    coll,
    name="pLIN.app",
    icon=None,  # .icns required for a real macOS dock icon; PNG isn't directly usable here, see build notes
    bundle_identifier="com.umcg.plin",
    info_plist={
        "NSHighResolutionCapable": "True",
        "CFBundleShortVersionString": PLIN_APP_VERSION,
        "CFBundleName": "pLIN",
        "NSRequiresAquaSystemAppearance": "False",
    },
)
