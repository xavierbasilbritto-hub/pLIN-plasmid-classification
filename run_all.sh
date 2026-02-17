#!/bin/bash
# ============================================================================
# pLIN — Complete Pipeline (macOS / Linux)
# Runs all steps: pLIN assignment → AMRFinderPlus → Integration
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
for dir in plasmid_sequences_for_training/*/fastas; do
    if [ -d "$dir" ]; then
        n=$(ls "$dir"/*.fasta 2>/dev/null | wc -l)
        FASTA_COUNT=$((FASTA_COUNT + n))
        inc=$(basename "$(dirname "$dir")")
        echo "  Found $n FASTA files in $inc"
    fi
done

if [ "$FASTA_COUNT" -eq 0 ]; then
    echo ""
    echo "  ERROR: No FASTA files found!"
    echo "  Place your plasmid FASTA files in:"
    echo "    plasmid_sequences_for_training/<IncGroup>/fastas/"
    exit 1
fi

echo ""
echo "  Total FASTA files: $FASTA_COUNT"
echo ""

# ── Step 1: Assign pLIN codes ─────────────────────────────────────────────
echo "============================================================"
echo "  STEP 1/3: Assigning pLIN codes"
echo "============================================================"
python3 assign_pLIN.py
echo ""

# ── Step 2: Run AMRFinderPlus ──────────────────────────────────────────────
echo "============================================================"
echo "  STEP 2/3: Running AMRFinderPlus"
echo "============================================================"

if command -v amrfinder &>/dev/null; then
    echo "  AMRFinderPlus detected. Running on all plasmids..."
    mkdir -p output/amrfinder
    for dir in plasmid_sequences_for_training/*/fastas; do
        if [ -d "$dir" ]; then
            for fasta in "$dir"/*.fasta; do
                [ -f "$fasta" ] || continue
                name=$(basename "$fasta" .fasta)
                out="output/amrfinder/${name}_amrfinder.tsv"
                if [ ! -f "$out" ]; then
                    amrfinder -n "$fasta" -o "$out" --plus 2>/dev/null || true
                fi
            done
        fi
    done
    echo "  AMRFinderPlus complete."
else
    echo "  WARNING: AMRFinderPlus not found. Skipping AMR detection."
    echo "  Install with: conda install -c bioconda -c conda-forge ncbi-amrfinderplus"
    echo ""
fi

# ── Step 3: Integrate pLIN + AMR ──────────────────────────────────────────
echo "============================================================"
echo "  STEP 3/3: Integrating pLIN with AMRFinderPlus"
echo "============================================================"

if [ -f "output/amrfinder/amrfinder_all_plasmids.tsv" ]; then
    python3 integrate_pLIN_AMR.py
else
    echo "  Skipping (AMRFinderPlus output not found)."
    echo "  Run AMRFinderPlus first, then: python3 integrate_pLIN_AMR.py"
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
echo ""
echo "  Launch the GUI:"
echo "    streamlit run plin_app.py"
echo ""
echo "============================================================"
