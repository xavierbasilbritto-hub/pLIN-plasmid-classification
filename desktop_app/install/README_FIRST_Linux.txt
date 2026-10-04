pLIN desktop app for Linux: installation and first use
=======================================================

pLIN gives every plasmid a permanent six-level code, from its backbone family (L1)
to its near-identical outbreak clone (L6). This folder contains the pLIN app for
Linux; no Python or other installation is needed.

What you need
-------------
- 64-bit Linux (x86_64) from 2024 or later, for example Ubuntu 24.04 or newer
  (the app needs glibc 2.39 or newer; on older systems install pLIN from the source code,
  see QUICK_START_Linux.md on GitHub)
- About 20 GB of free disk space (app about 1.5 GB, database 2 GB, search index 12 to 16 GB)
- An internet connection the first time (to download the database)

Step 1. Unzip
-------------
In a terminal, in the folder where you downloaded the file:
    unzip pLIN-linux.zip -d pLIN-app
    cd pLIN-app
The folder contains:
    pLIN/              the app
    sample_data/       example files to try the app (Swiss VIM-1 outbreak)
    README_FIRST.txt   this file

Step 2. Start pLIN
------------------
    chmod +x pLIN/pLIN
    ./pLIN/pLIN
After about 30 seconds pLIN opens in your web browser. Keep the terminal open while you
use pLIN. If the browser does not open, open the address shown in the terminal
("Starting local server at http://localhost:8501") yourself.

Step 3. Download the database (once)
------------------------------------
1. Keep "pLIN code scheme" on "v4.1 (recommended)".
2. Click "Download the pLIN v4.1 database (1.6 GB)". This takes a few minutes.
   Every file is checked against its published checksum. The database is saved in
   ~/.plin/plin_v41 and is used from then on.

Step 4. Try the example
-----------------------
1. Click "Browse files" under "Upload plasmid FASTA files" and select all 8 files in
   sample_data/swiss_vim1_outbreak (NARACHVIM11_plasmids.fasta ... NARACHVIM56_plasmids.fasta).
2. Click "Run Analysis".
   The first analysis builds a search index once (12 to 16 GB, 10 to 30 minutes).
   Later analyses take about a minute.
3. Open the Results tab and compare the pLIN codes with
   sample_data/swiss_vim1_outbreak/expected_pLIN_results.tsv.

Step 5. Analyse your own plasmids
---------------------------------
Upload FASTA files of complete plasmids or assembled plasmid contigs (.fasta, .fa, .fna).
Analyse all isolates of one investigation together. Export tables from the Export tab;
every export records the app and database version.

Quitting pLIN
-------------
Click "Quit pLIN" in the left sidebar, or press Ctrl+C in the terminal.

Good to know
------------
- The same sequences give the same codes on any computer. Codes of database plasmids,
  and every level that already exists in the database, never change. Parts of a code
  that are new to the database are marked "provisional" and are comparable only within
  one analysis.
- MMseqs2 and AMRFinderPlus are included in the Linux app.
- Problems: see the user manual (desktop_app/USER_MANUAL.md) at
  https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification
  or contact b.b.xavier@umcg.nl
