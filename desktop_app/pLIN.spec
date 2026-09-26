# -*- mode: python ; coding: utf-8 -*-
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
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

import os
from PyInstaller.utils.hooks import collect_data_files, copy_metadata, collect_submodules

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(SPEC)), ".."))

block_cipher = None

datas = []
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
    (os.path.join(REPO_ROOT, "data", "inc_classifier.npz"), "data"),
    (os.path.join(REPO_ROOT, "data", "inc_centroids.npz"), "data"),
    (os.path.join(REPO_ROOT, "output", "pLIN_assignments.tsv"), "output"),
    (os.path.join(REPO_ROOT, "assets"), "assets"),
    (os.path.join(REPO_ROOT, ".streamlit", "config.toml"), ".streamlit"),
]

hiddenimports = []
hiddenimports += collect_submodules("streamlit")
hiddenimports += collect_submodules("Bio")
hiddenimports += [
    "sklearn.utils._typedefs",
    "sklearn.neighbors._partition_nodes",
    "sklearn.utils._heap",
    "sklearn.utils._sorting",
    "sklearn.utils._vector_sentinel",
    "Bio.SeqIO",
    "Bio.Seq",
    "Bio.SeqRecord",
]
datas += collect_data_files("Bio")

a = Analysis(
    ["launcher.py"],
    pathex=[os.path.dirname(os.path.abspath(SPEC))],
    binaries=[],
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
    icon=None,  # .icns required for a real macOS dock icon; PNG isn't directly usable here — see build notes
    bundle_identifier="com.umcg.plin",
    info_plist={
        "NSHighResolutionCapable": "True",
        "CFBundleShortVersionString": "3.1.0",
        "CFBundleName": "pLIN",
        "NSRequiresAquaSystemAppearance": "False",
    },
)
