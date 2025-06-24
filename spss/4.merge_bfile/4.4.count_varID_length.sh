#!/bin/bash

FILELIST_DIR=/mnt/data1/lawless/spss/Samira.SPSS.data/filelist
PLINKOUT=/mnt/data1/lawless/spss/plink
OUTDIR=$PLINKOUT/count_varID_length

mkdir -p "$OUTDIR"

for CHR in {1..22}; do
  FILELIST="$FILELIST_DIR/chr${CHR}_files"
  while read -r line; do
    awk '{ print $2 " = " length($2) }' "$PLINKOUT/$line.bim" > "$OUTDIR/$line.count"
  done < "$FILELIST"
done

