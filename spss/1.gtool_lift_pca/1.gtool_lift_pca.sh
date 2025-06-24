#!/bin/bash

GTOOL=/home/lawless/tool/gtool
SAMPLE=/home/lawless/spss/SPSS.QC.sample.txt
DIR=/mnt/data1/lawless/spss/Samira.SPSS.data
PLINKOUT=/home/lawless/spss/plink
LOGDIR=$PLINKOUT/gtool_lift_pca_logs

mkdir -p "$LOGDIR"

for CHR in {1..22}; do
  FILELIST=$DIR/filelist/chr${CHR}_files
  LOGFILE=$LOGDIR/gtool_lift_pca.chr.${CHR}.log

  {
    while read gen; do
      $GTOOL -G \
        --g $DIR/$gen \
        --s $SAMPLE \
        --map $PLINKOUT/$gen.map \
        --ped $PLINKOUT/$gen.ped
    done < "$FILELIST"
  } 2>&1 | tee "$LOGFILE"
done

