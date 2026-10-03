# pLIN Desktop App: User Manual

pLIN (plasmid Lineage Identification Number) assigns every plasmid a
permanent, six-level code, from its backbone family (shared protein
families) to its near-identical outbreak clone (shared k-mers), so the same
plasmid gets the same code today and after the reference database has grown. This manual covers the **standalone desktop
app**: download it, run it, no separate Python install required.

If you're running pLIN from source (`streamlit run plin_app.py`) instead
of the standalone app, skip straight to [Using the app](#using-the-app):
everything from there on applies either way.

---

## 1. Download and install

Downloads are attached to each [GitHub Release](../../releases) as three
platform-specific zip files: `pLIN-macos.zip`, `pLIN-windows.zip`,
`pLIN-linux.zip`.

### macOS

1. Download `pLIN-macos.zip` and unzip it: you'll get `pLIN.app`.
2. Move `pLIN.app` to your Applications folder (or anywhere you like).
3. **First launch only**: macOS will refuse to open an app from an
   unidentified developer by default. Right-click (or Control-click)
   `pLIN.app` → **Open** → **Open** again in the dialog that appears.
   After this first time, double-clicking works normally.
4. A terminal window opens showing startup logs, and your default browser
   opens automatically to the app. **Leave the terminal window open**:
   closing it quits pLIN. To quit, close that window or press `Ctrl+C`
   inside it.

### Windows

1. Download `pLIN-windows.zip` and unzip it: you'll get a `pLIN` folder
   containing `pLIN.exe` and its supporting files. Keep the whole folder
   together; don't move `pLIN.exe` out on its own.
2. Double-click `pLIN.exe`.
3. **Windows Defender SmartScreen may show "Windows protected your PC"**
   the first time, since the app isn't code-signed. Click **More info**,
   then **Run anyway**.
4. A console window opens with startup logs, and your default browser
   opens automatically to the app. Leave the console window open while
   using pLIN; closing it quits the app.

### Linux

1. Download `pLIN-linux.zip` and unzip it: you'll get a `pLIN` folder.
2. Make the binary executable if needed: `chmod +x pLIN/pLIN`
3. Run it from a terminal: `./pLIN/pLIN`
4. Your default browser opens automatically to the app. Leave the
   terminal open; `Ctrl+C` there quits the app.

### If your browser doesn't open automatically

Look for a line in the terminal/console window like:

```
Starting local server at http://localhost:8501
```

Open that address manually in any browser (Chrome, Firefox, Safari,
Edge). pLIN runs entirely on your own machine, nothing is uploaded
anywhere unless you explicitly use a feature that says otherwise.

---

## 2. What's bundled, and what needs a separate install

The standalone app bundles the **core pipeline**: replicon (Inc/Rep group)
classification, pLIN code assignment, the cladogram viewer, within-group
sequence alignment, and outbreak/clone detection.

**pLIN v4.1 in the desktop app:**

1. **MMseqs2 is included.** The app bundles MMseqs2 (release 18-8cc5c, the version that built the
   database), which matches the proteins of new plasmids to the protein-family catalogue. On Windows,
   the first search sets up MMseqs2's helper tools and may ask once for administrator permission.
2. **The v4.1 database is downloaded once** (about 1.6 GB). With the code scheme on v4.1, the app
   shows a **Download the pLIN v4.1 database** button if the database is missing. It saves the files in
   `.plin/plin_v41` in your home folder (`~/.plin/plin_v41` on macOS/Linux,
   `%USERPROFILE%\.plin\plin_v41` on Windows) and checks every file against the SHA-256 checksums in
   `DATABASE_VERSION.json`. You can also download the files yourself from the GitHub
   [Releases](https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification/releases) page into
   that folder.
3. **A search index is built once.** The first time a plasmid with proteins new to the database is
   typed, the app builds a search index of about 16 GB in `~/.plin/v41_mmseqs_index` (10 to 30 minutes,
   depending on the computer).

Without them, the app tells you so and assigns the earlier v3 (legacy) codes instead.

To see exactly which database snapshot your copy has (plasmid count,
build date, unique pLIN codes), check the Export tab's caption after
running an analysis, or open `DATABASE_VERSION.json` inside the app's
installation folder. See [Will two colleagues get the same
code](#will-two-colleagues-get-the-same-code-for-the-same-plasmid) below
for why this matters, and the main repository README's "Database version
& updates" section for how to update to a newer database without
reinstalling the whole app.

A few **optional modules** call external bioinformatics tools that are
not bundled (they're large, platform-specific scientific packages
normally installed via [conda/bioconda](https://bioconda.github.io/)).
If a tool isn't installed, the app detects this automatically and either
hides that option or shows an install hint, nothing breaks, you simply
don't see that capability until the tool is present:

| Feature | Needs | Install |
|---|---|---|
| AMR gene screening | AMRFinderPlus | `conda install -c bioconda ncbi-amrfinderplus` |
| Mobility typing | MOB-suite | `conda install -c bioconda mob_suite` |
| True ANI validation | FastANI | `conda install -c bioconda fastani` |
| SNP sub-typing | minimap2 | `conda install -c bioconda minimap2` |
| CRISPR host inference | MinCED + BLAST+ | `conda install -c bioconda minced blast` |
| Chromosomal MLST typing | mlst | `conda install -c bioconda mlst` |
| Nucleotide Transformer (LLM) | transformers, torch | `pip install transformers torch` (not included in the standalone build; run from source if you need this) |
| Contrastive-encoder classifier (optional second opinion) | torch, plus a one-time training step | `pip install torch`, then run `python train_inc_encoder.py` once from a source checkout of the repository (this trains and saves the encoder artifact; it does not need to be repeated after that) |

If you install one of these tools after already launching pLIN, quit and
restart the app so it can detect the newly-installed tool.

### A note on the classifier choice

pLIN's default classifier is KNN (k-nearest-neighbour on raw tetranucleotide
composition, 91.1% cross-validated accuracy): fully interpretable, and what
the accompanying manuscript's validation is built on. When you select
"Auto-detect" for the Incompatibility Group, a **Classifier** option appears
letting you additionally try a contrastive-encoder classifier as a second
opinion when you have doubts about a specific call. Independent, leak-free
cross-validation shows the encoder modestly outperforms KNN overall (92.7%
accuracy, macro-F1 0.692 vs. KNN's 0.666), improving 19 of 28 replicon
groups and declining slightly on 8 (largest change: −0.016 F1, nothing
catastrophic). It never overrides the primary KNN result, it adds an
additional prediction alongside it. This option only appears if torch is
installed and the encoder has been trained (see the table above); it is
not required for normal use.

---

## 3. Using the app

### Step 1: Upload your plasmid sequences

On the **Overview** tab, upload one or more plasmid FASTA files
(`.fasta`, `.fa`, `.fna`) using the file uploader. You can upload several
plasmids at once. An optional metadata file (CSV/TSV) can be attached
alongside them if you want extra columns (collection date, source, etc.)
carried through to the results.

### Step 2: Choose optional analyses

Below the uploader, checkboxes appear for whichever optional tools the
app detected on your system (see the table above), for example, "Run
CRISPR host inference" or "Run MLST typing." Leave these unchecked if you
just want the core pLIN code assignment; check them if you have the
matching tool installed and want that extra analysis.

### Step 3: Run the analysis

Click **▶ Run Analysis**. The app then, automatically:

1. Predicts each plasmid's Inc/Rep replicon group (k-nearest-neighbour
   classifier: 91.1% cross-validated accuracy, macro-F1 0.67, 28 groups;
   accuracy is high mainly for common groups)
2. Assigns the six-level pLIN v4.1 code: proteins are matched to the
   protein-family catalogue and a k-mer sketch is taken. A plasmid with a near
   relative in the database copies that relative's backbone and lineage;
   otherwise it is placed in the backbone levels by the founder rule and starts
   a new lineage. A plasmid already in the database gets back exactly its
   published code; levels not yet in the database are flagged as new
3. Builds a distance tree of the uploaded plasmids for the dendrogram and
   heatmap views (this tree does not change the codes)
5. Screens for AMR genes with AMRFinderPlus, if installed and enabled

Click **🔄 Clear & Reset** at any point to start over with new files.

### Step 4: Read the results, tab by tab

- **Overview**: summary metrics (plasmid count, unique pLIN codes,
  AMR detections) and an input-sequence quality report flagging any
  sequence that's too short, has excessive ambiguous bases, or otherwise
  looks unreliable, so you know which results to double-check.
- **Results**: the full per-plasmid table: pLIN code, replicon type,
  confidence score, AMR gene list, and (if run) mobility/NT predictions.
  Sortable and exportable.
- **Cladogram**: the hierarchical clustering dendrogram for your
  uploaded plasmids, with the six pLIN threshold cuts overlaid. Also
  includes **within-group sequence alignment**: pick any two plasmids
  that share a pLIN code and see exactly how much of their sequence
  aligns and at what identity, rather than trusting the composition-based
  grouping blindly.
- **AMR Analysis**: gene-level detail from AMRFinderPlus: which genes,
  how many per plasmid, and prevalence charts, if AMR screening was run.
- **Epidemiology**: plasmid mobility prediction (conjugative vs.
  mobilisable vs. non-mobilisable) and the **outbreak/clone detection**
  module: plasmids in your dataset sharing an identical full six-level
  pLIN code are flagged as a likely shared-plasmid cluster, risk-rated
  HIGH/MODERATE/INFO based on shared AMR gene burden. This is the module
  that turns a permanent code into an actionable surveillance signal:
  see the accompanying manuscript's Discussion for why this matters.

  **If your isolates span multiple bacterial species or multiple Inc/Rep
  types**, check **Run MLST typing** in Step 2 and upload each isolate's
  assembled host genome (chromosome) alongside its plasmid FASTA, the
  app can auto-split a combined chromosome+plasmid assembly into its
  plasmid and chromosome contigs first if you upload it as one file (see
  "Contig Classification" in the Results tab). This enables
  **Pathogen-Plasmid Integration (MLST + pLIN)**, which combines
  chromosomal sequence type (ST) with pLIN code to distinguish two
  distinct outbreak signatures: **clonal spread** (same ST, same pLIN:
  one strain moving between patients) versus **horizontal plasmid
  transfer** (different STs or species, same pLIN, the plasmid itself
  spreading between different bacterial hosts, a real and clinically
  important outbreak pattern in its own right, independent of Inc/Rep
  typing since pLIN clusters by whole-plasmid composition rather than
  replicon family alone). A "Horizontal Transfer" count greater than
  zero here is itself a suspicion signal worth escalating: it means the
  same plasmid lineage was found in genuinely different bacterial
  backgrounds, which single-species or single-Inc-type surveillance
  would not catch. This exact pattern: a shared plasmid lineage
  disseminating across species and replicon-type boundaries: is what
  the bundled Swiss VIM-1 sample dataset demonstrates (see "Sample data"
  above), and has been validated against 13 independent, published
  multi-species outbreak studies.

  **Genome-plasmid linking note:** the app links each uploaded genome
  file to its plasmid by matching shared filename text (e.g.
  `sample01_genome.fasta` links to `sample01_plasmid.fasta`), or via a
  metadata CSV/TSV with a genome-name column and a `plasmid_id` column
  if your filenames don't share a common identifier. If MLST typing
  completes but no "Pathogen-Plasmid Integration" section appears
  afterwards, check for a warning about this, it means typing succeeded
  but no genome could be matched to a plasmid, so double-check your
  filenames share an identifiable token or provide a metadata file.
- **CRISPR Host**: plasmid-host association inference via CRISPR spacer
  matching, if MinCED/BLAST+ are installed and host genomes were
  uploaded.
- **DRAGNOME Buddy**: an interactive assistant for querying your
  results in plain language (requires a local Ollama installation;
  optional).
- **Export**: download your results as TSV tables, PDF/PNG figures, a
  bundled ZIP, or a JSON report.

### Query mode

If you upload plasmids without enough data to build a fresh reference
set, or explicitly choose to compare against the built-in database, pLIN
runs in **query mode**: it assigns your plasmid's code via
nearest-neighbour lookup against the bundled reference database rather
than clustering your plasmids from scratch. The Overview tab tells you
clearly when this is happening. The **pLIN status** column in your
results (see section 4 below) tells you whether each result matched an
existing database entry or was assigned a new code.

#### Will two colleagues get the same code for the same plasmid?

**If the plasmid already matches something in the reference database:
yes, always: as long as you're both on the same database version.** The
database file the app reads from is a fixed snapshot, it isn't modified
by running the app, and it doesn't change between a morning run and an
evening run, or between your computer and a colleague's. The same DNA
sequence always produces the same proteins and k-mer sketch, which always
find the same database relatives, which always yield the same code. Two people running the same already-catalogued plasmid, on
different computers, at different times, get identical results, provided
both installations report the same `database_version` (see
`DATABASE_VERSION.json`, or the version shown in the Export tab). The
**app version** (e.g. `v4.1.0`) and the **database version** (e.g.
`db-2026.10.03`) are tracked separately and can differ even between two
installs on the same app version: one of you may have downloaded a
newer standalone database update without upgrading the app itself, or
vice versa. Check both before assuming two runs are comparable.

**If the plasmid is genuinely new: not yet in anyone's reference
database: the two of you can get *different* numbers**, even for the
exact same sequence. A brand-new code isn't looked up; it's minted on
the spot as "one higher than the highest number this app instance
currently sees in its own copy of the database," entirely in memory, for
that run only. It is never written back to the database or shared with
anyone else's copy of the app. So:

- Running it twice on the *same* computer, same app install, same day,
  gives the same new number both times (nothing else changed in between).
- Running it on two different computers, or after either of you updates
  to a newer app release with an expanded database, can give two
  *different* new numbers for what is biologically the identical
  plasmid: because each app instance is counting up from its own
  locally-held database snapshot, independently, with no coordination
  between installs.

**Practical rule of thumb:** treat a matching, "Existing pLIN" result as
a permanent, shareable identifier you can compare across your team and
across time without hesitation. Treat a "New pLIN" result as
provisional and specific to that one run, useful for seeing how
divergent a plasmid is from anything currently known, but not yet a
stable label to compare between colleagues or sessions. If two of you
need to confirm whether your own novel outbreak strains are the same
lineage as each other, the reliable way is to compare the `nn_plasmid`,
`nn_distance`, and `new_plin_level` columns for each result (which
reference plasmid each one matched, how closely, and at what hierarchy
level it diverged) rather than comparing the raw new pLIN numbers
directly: or, simplest of all, upload both sequences together in the
same run, so they share one counter and get consistent numbers relative
to each other. A newly-discovered lineage only becomes a truly
permanent, shareable pLIN code once it's incorporated into a future
official database update.

#### Does switching classifiers (KNN vs. the optional encoder) change the pLIN code?

No. The pLIN code is computed from the plasmid's proteins and k-mers and
never looks at either classifier's output, so it is identical either way
for the same input sequence: only the reported Inc/Rep group label can
differ. See `pLIN_FAQ.txt` (question 2) at the root of the repository /
app installation folder for the full explanation.

---

## 4. Interpreting a pLIN code

A pLIN code has six dot-separated levels, from broadest to most specific:

```
L1.L2.L3.L4.L5.L6         (e.g. 169.178.183.208.209.2438, from the Swiss VIM-1 sample data)
│  │  │  │  │  └── L6: near-identical / outbreak clone (k-mer similarity >= 0.95 to its nearest relative)
│  │  │  │  └───── L5: lineage (k-mer similarity >= 0.80 to its nearest relative)
│  │  │  └──────── L4: backbone variant (k-mer containment >= 0.80)
│  │  └─────────── L3: shared backbone (k-mer containment >= 0.50)
│  └────────────── L2: backbone group (protein-family containment >= 0.60)
└───────────────── L1: backbone family (protein-family containment >= 0.40)
```

Two plasmids sharing all six levels are near-identical. Sharing the first
five means the same lineage. Sharing only the first two to four levels means
a related backbone with different gene content. The numbers are identifiers,
not distances. Existing codes never change when the database grows.

**pLIN status column.** Every query-mode result also reports whether its
code is "Existing pLIN" (the plasmid's full six-level code already
matches a plasmid already in the reference database) or "New pLIN" (it
required at least one freshly-minted digit, with `new_plin_level` naming
the first, coarsest level at which it diverged from its nearest
database match). This distinction matters for reproducibility across
computers and over time: see "Will two colleagues get the same code for
the same plasmid?" in section 3 above.

---

## 5. Troubleshooting

**"pLIN.app is damaged and can't be opened" (macOS)**, this is Gatekeeper
being overly cautious about an unsigned app, not actual file corruption.
Use the right-click → Open method in step 3 above. If it still refuses,
run this once in Terminal (adjust the path if you moved the app):
`xattr -cr /Applications/pLIN.app`

**Nothing happens when I double-click / the browser never opens**, check
the terminal/console window for an error message. The most common cause
is another program already using port 8501; pLIN should pick a different
free port automatically and print the actual URL it's using, look for
the `Starting local server at http://localhost:XXXX` line and open that
address manually.

**"Run analysis" is greyed out or nothing happens after upload**, make
sure your files are valid FASTA (`.fasta`/`.fa`/`.fna`) and not empty; the
Overview tab's input-quality report will tell you if a specific file was
rejected and why.

**An optional feature's checkbox is missing**: the corresponding tool
isn't installed or wasn't found. See the table in [section 2](#2-whats-bundled-and-what-needs-a-separate-install)
for the exact install command, then quit and restart pLIN.

**The app is very large / slow to start the first time**, this is
normal for a bundled scientific-Python application (numpy, pandas, scipy,
scikit-learn, matplotlib, and Streamlit's own frontend are all included).
Startup after the first launch is typically faster once your OS has
cached the files.

---

## 6. Getting help

- Source code, training data, and this manual:
  <https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification>
- To report a bug or request a feature, open a
  [GitHub Issue](../../issues) with your OS, pLIN version, and (if
  relevant) the exact error message from the terminal/console window.
- Citation is required if you use pLIN in published work, see
  `CITATION.cff` in the repository.
