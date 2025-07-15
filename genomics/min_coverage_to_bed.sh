#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# extract_high_cov_regions.sh
#
# Given an input BAM, extract regions with coverage >= MIN_COV and width > MIN_WIDTH
# Requires: bedtools, awk
#
# Usage:
#   ./extract_high_cov_regions.sh <input.bam> <output.bed>
# =============================================================================

# Minimum coverage threshold
MIN_COV=10

# Minimum region width (bp)
MIN_WIDTH=100

# Check arguments
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <input.bam> <output.bed>"
  exit 1
fi

IN_BAM="$1"
OUT_BED="$2"

# Verify input exists
if [[ ! -e "$IN_BAM" ]]; then
  echo "Error: input BAM '$IN_BAM' not found."
  exit 1
fi

# Verify output filename is non-empty
if [[ -z "$OUT_BED" ]]; then
  echo "Error: output filename is empty."
  exit 1
fi

# Extract high-coverage regions > MIN_WIDTH
bedtools genomecov -ibam "$IN_BAM" -bga \
  | awk -v cov="$MIN_COV" '$4 >= cov { print $1"\t"$2"\t"$3 }' \
  | bedtools sort -i - \
  | bedtools merge -i - \
  | awk -v w="$MIN_WIDTH" '$3 - $2 > w' \
  > "$OUT_BED"

echo "Wrote regions with coverage >= $MIN_COV and width > $MIN_WIDTH to '$OUT_BED'."
