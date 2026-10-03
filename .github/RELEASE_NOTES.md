## pLIN desktop app

pLIN gives every plasmid a permanent six-level code, from backbone family (L1) to near-identical outbreak clone (L6). This release runs **pLIN v4.1** codes on Windows, macOS and Linux.

### Download

| System | File |
|---|---|
| Windows 10/11 | `pLIN-windows.zip`: unzip and run `pLIN\pLIN.exe` |
| macOS | `pLIN-macos.zip`: unzip and open `pLIN.app` (first time: right-click, Open) |
| Linux | `pLIN-linux.zip`: unzip and run `pLIN/pLIN` |

The app opens in your web browser; keep its terminal window open while you work.

### First use

1. Keep the code scheme on **v4.1** and click **Download the pLIN v4.1 database (1.6 GB)**. Every file is checked against its published SHA-256 checksum and saved in `~/.plin/plin_v41` (Windows: `%USERPROFILE%\.plin\plin_v41`).
2. The first time a plasmid with proteins new to the database is typed, a search index is built once (about 12 to 16 GB of disk, a few minutes to half an hour).
3. To try the app, upload the 8 FASTA files in `sample_data/swiss_vim1_outbreak` (included) and compare with `expected_pLIN_results.tsv`.

MMseqs2 (18-8cc5c) is included. On Windows, the first search may ask once for administrator permission to set up MMseqs2's helper tools.

### Reproducibility

The same sequences give the same codes on any computer. Codes of database plasmids, and every level that exists in the database, never change. Parts of a code that are new to the database are marked as provisional (`provisional_from`) and are comparable within one analysis; analyse isolates of one investigation together. Every export records the app and database version.

Database release: `db-2026.10.03` (127,517 plasmids). Documentation: README and `desktop_app/USER_MANUAL.md`.
