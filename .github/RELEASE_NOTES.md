## pLIN desktop app

pLIN gives every plasmid a permanent six-level code, from backbone family (L1) to near-identical outbreak clone (L6). This release runs **pLIN v4.1** codes on Windows, macOS and Linux, with no Python or conda installation.

### Download

| System | File | Requirements |
|---|---|---|
| Windows | `pLIN-windows.zip` | Windows 10 or 11, 64-bit |
| macOS | `pLIN-macos.zip` | Mac with Apple silicon (M1 or later) |
| Linux | `pLIN-linux.zip` | 64-bit Linux from 2024 or later (e.g. Ubuntu 24.04) |

Each zip contains the app, **`README_FIRST.txt`** with step-by-step instructions for that system, and example data (`sample_data/swiss_vim1_outbreak`). Allow about 20 GB of free disk space.

### Quick steps

1. Unzip and start the app:
   - **Windows:** Extract All, open the `pLIN` folder, double-click `pLIN.exe` (at "Windows protected your PC": More info, Run anyway).
   - **macOS:** drag `pLIN.app` to Applications; the first time, right-click it and choose Open, then Open (or System Settings > Privacy & Security > Open Anyway).
   - **Linux:** `chmod +x pLIN/pLIN` and run `./pLIN/pLIN`.
2. pLIN opens in your web browser. Click **Download the pLIN v4.1 database (1.6 GB)** (once; every file is checked against its published checksum).
3. Try the example: upload the 8 files in `sample_data/swiss_vim1_outbreak`, click **Run Analysis**, and compare with `expected_pLIN_results.tsv`. The first analysis builds a search index once (12 to 16 GB, several minutes).
4. To stop pLIN, click **Quit pLIN** in the left sidebar.

MMseqs2 (18-8cc5c) is included. On Windows, the first MMseqs2 search may ask once for administrator permission to set up its helper tools.

### Reproducibility

The same sequences give the same codes on any computer. Codes of database plasmids, and every level that exists in the database, never change. Parts of a code that are new to the database are marked as provisional (`provisional_from`) and are comparable within one analysis, so analyse the isolates of one investigation together. Every export records the app and database version.

### Changes in 4.1.2

Licence is now the plain GNU General Public License v3.0 or later (the earlier mandatory-citation clause is replaced by a citation request in CITATION.cff). The typing method, thresholds and database are unchanged, so codes are identical to 4.1.1. Adds the scripts of the external validation and sensitivity analyses.

### Changes in the database (db-2026.10.05)

The replicon classifier gains the **IncL/M** group. Without it, pOXA-48-like plasmids were assigned to the nearest available type, usually IncFII; this affected 1,352 plasmids including the largest epidemic clone in the database. The classifier now covers 29 groups (8,404 training plasmids, 91.0% cross-validated accuracy, IncL/M F1 0.98) and the `inc_type` column of the release has been regenerated with it.

**pLIN codes are unchanged.** Every code, sketch, protein family and sequence hash is identical to db-2026.10.03, so codes issued from the previous release stay valid; only the replicon annotation changed.

Database release: `db-2026.10.05` (127,517 plasmids). Archived on Zenodo: https://doi.org/10.5281/zenodo.23126056 (concept DOI, resolves to the latest version). Full instructions: `desktop_app/USER_MANUAL.md` and the README.
