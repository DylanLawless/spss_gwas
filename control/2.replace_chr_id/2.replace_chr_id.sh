#!/bin/bash

FILELIST_DIR=/mnt/data1/lawless/spss/gwas/filelist
PLINKOUT=/mnt/data1/lawless/spss/gwas/shcs

for CHR in {1..22}; do
  while read -r line; do
    awk -v chr="$CHR" '{$1 = chr; print}' "$PLINKOUT/$line.map" > "$PLINKOUT/$line.chr.map"
  done < "$FILELIST_DIR/chr${CHR}_files"
done

