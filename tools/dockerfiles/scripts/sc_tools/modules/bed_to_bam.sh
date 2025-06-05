#!/bin/bash
set -euo pipefail
BED_FILE="$1"
BAM_FILE="$2"
CHROM_SIZES="$3"
THREADS="$4"
TEMP_BAM="${BAM_FILE%.bam}_unsorted.bam"
bedtools bedtobam -i "$BED_FILE" -g "$CHROM_SIZES" > "$TEMP_BAM"
samtools sort -@ "$THREADS" -o "$BAM_FILE" "$TEMP_BAM"
samtools index "$BAM_FILE"
rm -f "$TEMP_BAM"