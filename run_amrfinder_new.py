#!/usr/bin/env python3
"""Run AMRFinderPlus on the 17 new Inc groups and append to existing results."""
import subprocess, os, sys, glob

AMRFINDER = "/Users/basilxavier/miniconda3/envs/pLIN_analysis/bin/amrfinder"
BASE = "/Users/basilxavier/Desktop/PLASMID_TOOL/plasmid_sequences_for_training"
EXISTING_FILE = "/Users/basilxavier/Desktop/PLASMID_TOOL/output/amrfinder/amrfinder_all_plasmids.tsv"

EXISTING_GROUPS = {"IncFII", "IncN", "IncX1"}
ALL_GROUPS = sorted([d for d in os.listdir(BASE) if os.path.isdir(os.path.join(BASE, d))])
NEW_GROUPS = [g for g in ALL_GROUPS if g not in EXISTING_GROUPS]

print(f"Processing {len(NEW_GROUPS)} new Inc groups...", flush=True)

all_new_rows = []
total_processed = 0
total_hits = 0
errors = 0

for inc in NEW_GROUPS:
    fasta_dir = os.path.join(BASE, inc, "fastas")
    fastas = sorted(glob.glob(os.path.join(fasta_dir, "*.fasta")))
    print(f"\n{inc}: {len(fastas)} sequences", flush=True)

    for i, fasta_path in enumerate(fastas):
        fname = os.path.basename(fasta_path).replace(".fasta", "")
        total_processed += 1

        try:
            result = subprocess.run(
                [AMRFINDER, "-n", fasta_path, "--plus",
                 "--threads", "4", "--name", fname],
                capture_output=True, text=True, timeout=120
            )
            if result.returncode != 0 and "ERROR" in result.stderr:
                print(f"  ERROR {fname}: {result.stderr.strip().split(chr(10))[-1]}", flush=True)
                errors += 1
                continue

            lines = result.stdout.strip().split("\n")
            data_lines = [l for l in lines[1:] if l.strip()]

            for line in data_lines:
                all_new_rows.append(f"{fname}\t{inc}\t{line}")
                total_hits += 1

        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT: {fname}", flush=True)
            errors += 1
        except Exception as e:
            print(f"  ERROR {fname}: {e}", flush=True)
            errors += 1

        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{len(fastas)} done ({total_hits} hits so far)", flush=True)

    print(f"  Done. Running total: {total_hits} hits from {total_processed} sequences", flush=True)

print(f"\nTotal: {total_hits} new AMR hits from {total_processed} sequences ({errors} errors)", flush=True)
if all_new_rows:
    with open(EXISTING_FILE, "a") as f:
        for row in all_new_rows:
            f.write(row + "\n")
    print(f"Appended {len(all_new_rows)} rows to {EXISTING_FILE}", flush=True)
else:
    print("No AMR hits found in new groups", flush=True)

print(f"\nDone! Processed {total_processed} sequences across {len(NEW_GROUPS)} Inc groups", flush=True)
