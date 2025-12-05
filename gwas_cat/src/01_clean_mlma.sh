#!/usr/bin/env bash

INFILE="$1"
OUTFILE="${INFILE%.txt}_snps_only.tsv"
OUTFILE_NOSNP="${INFILE%.txt}_snps_only_noSNP.tsv"

if [ -z "$INFILE" ]; then
  echo "Usage: ./clean_mlma.sh ../data/cohort.filt.group.case_control.loco.mlma.printed.sig.txt"
  exit 1
fi

echo "Cleaning GCTA MLMA file: $INFILE"
echo "Temporary SNP cleaned file: $OUTFILE"
echo "Final file without SNP column: $OUTFILE_NOSNP"

A1_COL=4
A2_COL=5

# Step 1: keep only SNPs (A,C,G,T alleles)
awk -v a1="$A1_COL" -v a2="$A2_COL" '
BEGIN { FS=OFS="\t" }
NR==1 { print; next }
($a1 ~ /^[ACGT]$/ && $a2 ~ /^[ACGT]$/) { print }
' "$INFILE" > "$OUTFILE"

# Step 2: remove SNP column entirely (column 2)
cut -f1,3- "$OUTFILE" > "$OUTFILE_NOSNP"

echo "Done."
echo "Remaining variant count:"
wc -l "$OUTFILE_NOSNP"

