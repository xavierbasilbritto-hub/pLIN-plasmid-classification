#!/bin/bash
# Run AMRFinderPlus on all plasmid sequences and produce a combined output
# Uses the pLIN_analysis conda environment's AMRFinderPlus v4.2.5

set -e

AMRFINDER="/Users/basilxavier/miniconda3/envs/pLIN_analysis/bin/amrfinder"
DB_DIR="/Users/basilxavier/miniconda3/envs/pLIN_analysis/share/amrfinderplus/data/2026-01-21.1"
BASE_DIR="/Users/basilxavier/Desktop/PLASMID_TOOL"
OUTPUT_DIR="$BASE_DIR/output/amrfinder"
TRAINING_DIR="$BASE_DIR/plasmid_sequences_for_training"

mkdir -p "$OUTPUT_DIR"

# Combined output file with header
COMBINED="$OUTPUT_DIR/amrfinder_all_plasmids.tsv"
echo -e "source_file\tinc_type\tProtein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description" > "$COMBINED"

total=0
with_hits=0
processed=0

for INC_TYPE in IncFII IncN IncX1; do
    FASTA_DIR="$TRAINING_DIR/$INC_TYPE/fastas"
    echo "Processing $INC_TYPE ..."
    count=0

    for FASTA in "$FASTA_DIR"/*.fasta; do
        BASENAME=$(basename "$FASTA" .fasta)

        # Run AMRFinderPlus (nucleotide mode with --plus for virulence/stress genes)
        RESULT=$($AMRFINDER -n "$FASTA" -d "$DB_DIR" --plus 2>/dev/null | tail -n +2)

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
