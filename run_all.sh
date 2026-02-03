#!/bin/bash
# ============================================================================
# pLIN — Complete Pipeline (macOS / Linux)
# Runs all steps: pLIN assignment → AMRFinderPlus → Integration → Figures
# ============================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "============================================================"
echo "  pLIN — Complete Analysis Pipeline"
echo "============================================================"
echo ""

# ── Activate virtual environment ───────────────────────────────────────────
if [ -d ".venv" ]; then
    source .venv/bin/activate
    echo "  Using Python: $(python3 --version)"
else
    echo "  WARNING: .venv not found. Run setup.sh first."
    echo "  Continuing with system Python ..."
fi

# ── Check for FASTA files ─────────────────────────────────────────────────
FASTA_COUNT=0
for inc in IncFII IncN IncX1; do
    dir="plasmid_sequences_for_training/$inc/fastas"
    if [ -d "$dir" ]; then
        n=$(ls "$dir"/*.fasta 2>/dev/null | wc -l)
        FASTA_COUNT=$((FASTA_COUNT + n))
        echo "  Found $n FASTA files in $inc"
    fi
done

if [ "$FASTA_COUNT" -eq 0 ]; then
    echo ""
    echo "  ERROR: No FASTA files found!"
    echo "  Place your plasmid FASTA files in:"
    echo "    plasmid_sequences_for_training/IncFII/fastas/"
    echo "    plasmid_sequences_for_training/IncN/fastas/"
    echo "    plasmid_sequences_for_training/IncX1/fastas/"
    exit 1
fi

echo ""
echo "  Total FASTA files: $FASTA_COUNT"
echo ""

# ── Step 1: Assign pLIN codes ─────────────────────────────────────────────
echo "============================================================"
echo "  STEP 1/4: Assigning pLIN codes"
echo "============================================================"
python3 assign_pLIN.py
echo ""

# ── Step 2: Run AMRFinderPlus ──────────────────────────────────────────────
echo "============================================================"
echo "  STEP 2/4: Running AMRFinderPlus"
echo "============================================================"

if command -v amrfinder &>/dev/null || [ -n "$AMRFINDER_PATH" ]; then
    bash run_amrfinder_all.sh
else
    # Try to find it in conda envs
    FOUND_AMR=0
    for candidate in \
        "$HOME/miniconda3/envs/*/bin/amrfinder" \
        "$HOME/miniforge3/envs/*/bin/amrfinder" \
        "$HOME/anaconda3/envs/*/bin/amrfinder"; do
        for f in $candidate; do
            if [ -x "$f" ]; then
                FOUND_AMR=1
                bash run_amrfinder_all.sh
                break 2
            fi
        done
    done

    if [ "$FOUND_AMR" -eq 0 ]; then
        echo "  WARNING: AMRFinderPlus not found. Skipping AMR detection."
        echo "  Install with: conda install -c bioconda ncbi-amrfinderplus"
        echo "  Then re-run: bash run_amrfinder_all.sh"
        echo ""
    fi
fi

# ── Step 3: Integrate pLIN + AMR ──────────────────────────────────────────
echo "============================================================"
echo "  STEP 3/4: Integrating pLIN with AMRFinderPlus"
echo "============================================================"

if [ -f "output/amrfinder/amrfinder_all_plasmids.tsv" ]; then
    python3 integrate_pLIN_AMR.py
else
    echo "  Skipping (AMRFinderPlus output not found)."
    echo "  Run AMRFinderPlus first, then: python3 integrate_pLIN_AMR.py"
fi
echo ""

# ── Step 4: Generate figures ───────────────────────────────────────────────
echo "============================================================"
echo "  STEP 4/4: Generating publication figures"
echo "============================================================"

if [ -f "output/integrated/pLIN_AMR_integrated.tsv" ]; then
    python3 generate_figures.py
else
    echo "  Skipping (integrated data not found)."
    echo "  Complete steps 2-3 first, then: python3 generate_figures.py"
fi

echo ""
echo "============================================================"
echo "  Pipeline Complete!"
echo "============================================================"
echo ""
echo "  Output files:"
echo "    output/pLIN_assignments.tsv           — pLIN codes"
echo "    output/amrfinder/                     — AMR detections"
echo "    output/integrated/                    — Combined tables"
echo "    output/figures/                       — Publication figures"
echo ""
echo "============================================================"
