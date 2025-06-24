#!/bin/bash

FILELIST_DIR=/mnt/data1/lawless/spss/gwas/filelist
MERGE_LIST=/mnt/data1/lawless/spss/gwas/shcs/4.merge_bfile/allfiles.txt
mkdir -p "$(dirname "$MERGE_LIST")"
> "$MERGE_LIST"

for CHR in {1..22}; do
  INPUT="$FILELIST_DIR/chr${CHR}_files.txt"
  OUTPUT="$FILELIST_DIR/chr${CHR}_files.merge.txt"

  awk '{print $1".bed", $1".bim", $1".fam"}' "$INPUT" > "$OUTPUT"

  # skip chr1 in final list (used as base in --bfile)
  if [ "$CHR" -ne 1 ]; then
    cat "$OUTPUT" >> "$MERGE_LIST"
  fi
done

