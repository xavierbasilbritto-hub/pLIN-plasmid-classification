#!/bin/bash
# ============================================================================
# Run AMRFinderPlus on all plasmid sequences and produce a combined output
# Portable version — auto-detects AMRFinderPlus binary and database
# ============================================================================
set -e

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTPUT_DIR="$BASE_DIR/output/amrfinder"
TRAINING_DIR="$BASE_DIR/plasmid_sequences_for_training"

# ── Auto-detect AMRFinderPlus ──────────────────────────────────────────────
# Check common locations for the amrfinder binary
if [ -n "$AMRFINDER_PATH" ]; then
    AMRFINDER="$AMRFINDER_PATH"
elif command -v amrfinder &>/dev/null; then
    AMRFINDER="$(command -v amrfinder)"
elif [ -f "$CONDA_PREFIX/bin/amrfinder" ]; then
    AMRFINDER="$CONDA_PREFIX/bin/amrfinder"
else
    # Search common conda env locations
    for candidate in \
        "$HOME/miniconda3/envs/*/bin/amrfinder" \
        "$HOME/miniforge3/envs/*/bin/amrfinder" \
        "$HOME/anaconda3/envs/*/bin/amrfinder" \
        "$HOME/mambaforge/envs/*/bin/amrfinder"; do
        for f in $candidate; do
            if [ -x "$f" ]; then
                AMRFINDER="$f"
                break 2
            fi
        done
    done
fi

if [ -z "$AMRFINDER" ] || [ ! -x "$AMRFINDER" ]; then
    echo "ERROR: AMRFinderPlus not found."
    echo ""
    echo "Install it with:"
    echo "  conda install -c bioconda -c conda-forge ncbi-amrfinderplus"
    echo "  amrfinder -u    # update database"
    echo ""
    echo "Or set AMRFINDER_PATH environment variable:"
    echo "  export AMRFINDER_PATH=/path/to/amrfinder"
    exit 1
fi

echo "Using AMRFinderPlus: $AMRFINDER"
echo "Version: $($AMRFINDER --version 2>&1 | head -1)"

# ── Auto-detect database ──────────────────────────────────────────────────
if [ -n "$AMRFINDER_DB" ]; then
    DB_DIR="$AMRFINDER_DB"
else
    # Try to find the database in the same prefix as the binary
    AMRFINDER_PREFIX="$(dirname "$(dirname "$AMRFINDER")")"
    DB_BASE="$AMRFINDER_PREFIX/share/amrfinderplus/data"
    if [ -d "$DB_BASE" ]; then
        # Use the latest database version
        DB_DIR="$(ls -1d "$DB_BASE"/20* 2>/dev/null | sort -V | tail -1)"
    fi
fi

if [ -z "$DB_DIR" ] || [ ! -d "$DB_DIR" ]; then
    echo "WARNING: AMRFinderPlus database not found. Running without -d flag."
    echo "  If this fails, update the database with: amrfinder -u"
    DB_FLAG=""
else
    echo "Using database: $DB_DIR"
    DB_FLAG="-d $DB_DIR"
fi

# ── Check input data ──────────────────────────────────────────────────────
if [ ! -d "$TRAINING_DIR" ]; then
    echo "ERROR: Training directory not found: $TRAINING_DIR"
    echo "  Place FASTA files in plasmid_sequences_for_training/<IncType>/fastas/"
    exit 1
fi

echo ""
echo "========================================"
echo "AMRFinderPlus — Processing all plasmids"
echo "========================================"

mkdir -p "$OUTPUT_DIR"

# Combined output file with header
COMBINED="$OUTPUT_DIR/amrfinder_all_plasmids.tsv"
echo -e "source_file\tinc_type\tProtein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description" > "$COMBINED"

total=0
with_hits=0
processed=0

for INC_TYPE in IncFII IncN IncX1; do
    FASTA_DIR="$TRAINING_DIR/$INC_TYPE/fastas"

    if [ ! -d "$FASTA_DIR" ]; then
        echo "  Skipping $INC_TYPE (directory not found: $FASTA_DIR)"
        continue
    fi

    # Check if there are FASTA files
    shopt -s nullglob
    fasta_files=("$FASTA_DIR"/*.fasta)
    shopt -u nullglob

    if [ ${#fasta_files[@]} -eq 0 ]; then
        echo "  Skipping $INC_TYPE (no .fasta files found)"
        continue
    fi

    echo "Processing $INC_TYPE (${#fasta_files[@]} files) ..."
    count=0

    for FASTA in "${fasta_files[@]}"; do
        BASENAME=$(basename "$FASTA" .fasta)

        # Run AMRFinderPlus (nucleotide mode with --plus for virulence/stress genes)
        RESULT=$($AMRFINDER -n "$FASTA" $DB_FLAG --plus 2>/dev/null | tail -n +2) || true

        if [ -n "$RESULT" ]; then
            # Prepend source file and inc_type to each line
            echo "$RESULT" | while IFS= read -r line; do
                echo -e "${BASENAME}\t${INC_TYPE}\t${line}"
            done >> "$COMBINED"
            with_hits=$((with_hits + 1))
        fi

        total=$((total + 1))
        count=$((count + 1))
        processed=$((processed + 1))

        if [ $((processed % 500)) -eq 0 ]; then
            echo "  Processed $processed plasmids so far ($with_hits with hits) ..."
        fi
    done

    echo "  $INC_TYPE: $count plasmids processed"
done

echo ""
echo "========================================"
echo "AMRFinderPlus Complete"
echo "========================================"
echo "Total plasmids processed: $total"
echo "Plasmids with AMR/virulence hits: $with_hits"
echo "Output file: $COMBINED"
echo "========================================"
