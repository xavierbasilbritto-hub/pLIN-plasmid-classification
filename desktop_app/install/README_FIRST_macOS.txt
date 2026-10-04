pLIN desktop app for macOS: installation and first use
=======================================================

pLIN gives every plasmid a permanent six-level code, from its backbone family (L1)
to its near-identical outbreak clone (L6). This folder contains the pLIN app for
macOS; no Python or other installation is needed.

What you need
-------------
- A Mac with Apple silicon (M1, M2, M3, M4 or later). For Intel Macs, install pLIN from
  the source code instead (QUICK_START_macOS.md on GitHub).
- About 20 GB of free disk space (app about 1.5 GB, database 2 GB, search index 12 to 16 GB)
- An internet connection the first time (to download the database)

Step 1. Unzip and move the app
------------------------------
1. Double-click pLIN-macos.zip. A folder appears containing:
     pLIN.app           the app
     sample_data        example files to try the app (Swiss VIM-1 outbreak)
     README_FIRST.txt   this file
2. Drag pLIN.app into your Applications folder.

Step 2. Open pLIN the first time
--------------------------------
The app is not signed with an Apple developer certificate, so macOS asks once:
1. In Applications, right-click (or Control-click) pLIN.app and choose "Open",
   then click "Open" again.
2. If macOS only offers "Move to Trash" or "Done": click "Done", open
   System Settings > Privacy & Security, scroll down to the message about pLIN and
   click "Open Anyway", then confirm.
   (Alternative: in Terminal run   xattr -cr /Applications/pLIN.app   and open the app again.)
After this first time, double-clicking pLIN.app works normally.

Step 3. Use pLIN in your browser
--------------------------------
pLIN has no window of its own: after about 30 seconds it opens in your web browser
(Safari, Chrome or Firefox) at an address such as http://localhost:8501.
If no browser tab appears, open http://localhost:8501 yourself.

Step 4. Download the database (once)
------------------------------------
1. Keep "pLIN code scheme" on "v4.1 (recommended)".
2. Click "Download the pLIN v4.1 database (1.6 GB)". This takes a few minutes.
   Every file is checked against its published checksum. The database is saved in
   ~/.plin/plin_v41 (in your home folder) and is used from then on.

Step 5. Try the example
-----------------------
1. Click "Browse files" under "Upload plasmid FASTA files" and select all 8 files in
   sample_data/swiss_vim1_outbreak (NARACHVIM11_plasmids.fasta ... NARACHVIM56_plasmids.fasta).
2. Click "Run Analysis".
   The first analysis builds a search index once (12 to 16 GB, 3 to 30 minutes).
   Later analyses take about a minute.
3. Open the Results tab and compare the pLIN codes with
   sample_data/swiss_vim1_outbreak/expected_pLIN_results.tsv (open it in Numbers or Excel).

Step 6. Analyse your own plasmids
---------------------------------
Upload FASTA files of complete plasmids or assembled plasmid contigs (.fasta, .fa, .fna).
Analyse all isolates of one investigation together. Export tables from the Export tab;
every export records the app and database version.

Quitting pLIN
-------------
Click "Quit pLIN" in the left sidebar of the app. Closing the browser tab alone does not
stop pLIN. (If needed: open Activity Monitor, select pLIN and click Quit.)

Good to know
------------
- The same sequences give the same codes on any computer. Codes of database plasmids,
  and every level that already exists in the database, never change. Parts of a code
  that are new to the database are marked "provisional" and are comparable only within
  one analysis.
- MMseqs2 and AMRFinderPlus are included in the macOS app.
- Problems: see the user manual (desktop_app/USER_MANUAL.md) at
  https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification
  or contact b.b.xavier@umcg.nl
