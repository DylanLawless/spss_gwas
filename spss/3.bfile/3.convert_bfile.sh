#!/bin/bash

FILELIST_DIR=/mnt/data1/lawless/spss/Samira.SPSS.data/filelist
PLINKOUT=/mnt/data1/lawless/spss/plink
PLINK=/home/lawless/tool/plink

# make-bed for each file
for CHR in {1..22}; do
  while read -r line; do
    $PLINK --file "$PLINKOUT/$line" --make-bed --out "$PLINKOUT/$line"
    $PLINK --file "$PLINKOUT/$line" --recode --out "$PLINKOUT/${line}.recode"
  done < "$FILELIST_DIR/chr${CHR}_files"
done

