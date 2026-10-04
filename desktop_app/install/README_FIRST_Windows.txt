pLIN desktop app for Windows: installation and first use
=========================================================

pLIN gives every plasmid a permanent six-level code, from its backbone family (L1)
to its near-identical outbreak clone (L6). This folder contains the pLIN app for
Windows; no Python or other installation is needed.

What you need
-------------
- Windows 10 or 11 (64-bit)
- About 20 GB of free disk space (app 1 GB, database 2 GB, search index 12 to 16 GB)
- An internet connection the first time (to download the database)

Step 1. Unzip
-------------
1. Right-click pLIN-windows.zip and choose "Extract All...", then "Extract".
   Do not run the app from inside the zip file.
2. Open the extracted folder. It contains:
     pLIN\              the app (keep this folder together; do not move pLIN.exe out of it)
     sample_data\       example files to try the app (Swiss VIM-1 outbreak)
     README_FIRST.txt   this file

Step 2. Start pLIN
------------------
1. Open the pLIN folder and double-click pLIN.exe.
2. If Windows shows "Windows protected your PC", click "More info", then "Run anyway".
   (The app is not code-signed; this appears only the first time.)
3. A black console window opens and, after about 30 seconds, pLIN opens in your web
   browser. Keep the console window open while you use pLIN.
   If the browser does not open, look in the console window for a line such as
   "Starting local server at http://localhost:8501" and open that address in Chrome,
   Edge or Firefox.

Step 3. Download the database (once)
------------------------------------
1. In pLIN, keep "pLIN code scheme" on "v4.1 (recommended)".
2. Click "Download the pLIN v4.1 database (1.6 GB)". This takes a few minutes.
   Every file is checked against its published checksum. The database is saved in
   %USERPROFILE%\.plin\plin_v41 and is used from then on.

Step 4. Try the example
-----------------------
1. Click "Browse files" under "Upload plasmid FASTA files" and select all 8 files in
   sample_data\swiss_vim1_outbreak (NARACHVIM11_plasmids.fasta ... NARACHVIM56_plasmids.fasta).
2. Click "Run Analysis".
   The first analysis builds a search index once (12 to 16 GB, 10 to 30 minutes).
   The first time MMseqs2 runs, Windows may ask for administrator permission to set up
   its helper tools: click "Yes". Later analyses take about a minute.
3. Open the Results tab and compare the pLIN codes with
   sample_data\swiss_vim1_outbreak\expected_pLIN_results.tsv (open it in Excel).

Step 5. Analyse your own plasmids
---------------------------------
Upload FASTA files of complete plasmids or assembled plasmid contigs (.fasta, .fa, .fna).
Analyse all isolates of one investigation together. Export tables from the Export tab;
every export records the app and database version.

Quitting pLIN
-------------
Click "Quit pLIN" in the left sidebar, or close the black console window.

Good to know
------------
- The same sequences give the same codes on any computer. Codes of database plasmids,
  and every level that already exists in the database, never change. Parts of a code
  that are new to the database are marked "provisional" and are comparable only within
  one analysis.
- AMRFinderPlus (resistance genes) is optional and is not included for Windows. pLIN
  uses it if it is installed in WSL (Windows Subsystem for Linux). pLIN codes do not
  need it.
- Problems: see the user manual (desktop_app/USER_MANUAL.md) at
  https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification
  or contact b.b.xavier@umcg.nl
