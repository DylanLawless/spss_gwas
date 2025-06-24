#!/bin/bash

PLINK=/home/lawless/tool/plink
PLINKDIR=/mnt/data1/lawless/spss/plink
MERGEDIR=$PLINKDIR/replaced.ID.bim
MERGE_LIST=$PLINKDIR/4.merge_bfile/allfiles.txt
MYFILE1=SPSS.QC_chr1.pos13417-5012355.impute2

$PLINK \
  --bfile "$MERGEDIR/$MYFILE1" \
  --merge-list "$MERGE_LIST" \
  --make-bed \
  --out "$MERGEDIR/SPSS.QC.impute.allfiles" \
  --threads 26 \
  --memory 102400

