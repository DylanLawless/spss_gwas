#!/bin/bash

PLINKOUT=/mnt/data1/lawless/spss/plink
COUNTDIR=$PLINKOUT/count_varID_length

mkdir -p "$COUNTDIR/over79"

for CHR in {1..22}; do
  for file in "$COUNTDIR"/SPSS.QC_chr${CHR}.*.count; do
    awk -F '=' '$2 > 79 { print $0 }' "$file" > "${file}.varID_over79"
  done
done

