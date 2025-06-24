#!/bin/bash

PLINKOUT=/mnt/data1/lawless/spss/plink
COUNTDIR=$PLINKOUT/count_varID_length
REPLACED_DIR=$PLINKOUT/replaced.ID.bim

mkdir -p "$REPLACED_DIR"

for CHR in {1..22}; do
  for file in "$PLINKOUT"/SPSS.QC_chr${CHR}.*.impute2.bim; do
    sed -f "$COUNTDIR/All.varID_over79.sed.list" "$file" > "$file.replaced.ID"
    mv "$file.replaced.ID" "$REPLACED_DIR/$(basename "$file")"
  done
done

# Optional: move original bim files if needed
# mkdir -p "$PLINKOUT/original_bim"
# mv "$PLINKOUT"/SPSS.QC_chr*.impute2.bim "$PLINKOUT/original_bim/"

# Optional: move accompanying bed/fam files (assumes same basename)
# mv "$PLINKOUT"/SPSS.QC_chr*.impute2.{bed,fam} "$REPLACED_DIR/"

# Optional: strip .replaced.ID if you previously used such suffixes
# cd "$REPLACED_DIR"
# for file in *.bim; do mv "$file" "$(echo "$file" | sed 's/.replaced.ID//')"; done

