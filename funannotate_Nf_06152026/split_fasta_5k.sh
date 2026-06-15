#!/bin/bash

# Usage: ./split_fasta.sh input.fasta

if [ $# -ne 1 ]; then
    echo "Usage: $0 input.fasta"
    exit 1
fi

INPUT=$1
PREFIX=$(basename "$INPUT" .fasta)
COUNT=0
FILECOUNT=1

OUTFILE="${PREFIX}_${FILECOUNT}.fasta"
> "$OUTFILE"

while read -r LINE; do
    if [[ $LINE == ">"* ]]; then
        # New sequence header
        ((COUNT++))
        if (( COUNT > 5000 )); then
            # Start new file
            ((FILECOUNT++))
            OUTFILE="${PREFIX}_${FILECOUNT}.fasta"
            > "$OUTFILE"
            COUNT=1
        fi
    fi
    echo "$LINE" >> "$OUTFILE"
done < "$INPUT"

echo "Done. Created $FILECOUNT files."
