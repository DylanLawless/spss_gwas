#!/bin/bash

PLINKOUT=/mnt/data1/lawless/spss/gwas/shcs
COUNTDIR=$PLINKOUT/count_varID_length
REPLACED_DIR=$PLINKOUT/replaced.ID.bim

mkdir -p "$REPLACED_DIR"

for CHR in {1..22}; do
  for file in "$PLINKOUT"/SHCS.QC_chr${CHR}.*.impute2.bim; do
    sed -f "$COUNTDIR/All.varID_over79.sed.list" "$file" > "$REPLACED_DIR/$(basename "$file")"
  done
done

# Optional: move accompanying .bed and .fam files
for EXT in bed fam; do
  find "$PLINKOUT" -name "SHCS.QC_chr*.impute2.$EXT" -exec mv {} "$REPLACED_DIR/" \;
done

# Optional: move original .bim to backup dir
mkdir -p "$PLINKOUT/original_bim"
find "$PLINKOUT" -name "SHCS.QC_chr*.impute2.bim" -exec mv {} "$PLINKOUT/original_bim/" \;

