#!/bin/bash

GTOOL=/home/lawless/tool/gtool
SAMPLE=/mnt/data1/lawless/spss/gwas/SHCS.QC.sample.txt
DIR=/mnt/data1/lawless/spss/gwas/Samira.SHCS.data
PLINKOUT=/mnt/data1/lawless/spss/gwas/plink_hiv
LOGDIR=$PLINKOUT/gtool_lift_pca_logs

mkdir -p "$LOGDIR"

for CHR in {1..22}; do
  FILELIST=/mnt/data1/lawless/spss/gwas/filelist/chr${CHR}_files
  LOGFILE=$LOGDIR/gtool_lift_pca.chr.${CHR}.log

  {
    while read -r gen; do
      $GTOOL -G \
        --g $DIR/$gen \
        --s $SAMPLE \
        --map $PLINKOUT/$gen.map \
        --ped $PLINKOUT/$gen.ped
    done < "$FILELIST"
  } 2>&1 | tee "$LOGFILE"
done

