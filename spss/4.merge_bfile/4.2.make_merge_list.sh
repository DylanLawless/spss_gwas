#!/bin/bash

PLINKDIR=/mnt/data1/lawless/spss/plink
MERGE_LIST=$PLINKDIR/allfiles.txt

> "$MERGE_LIST"

for CHR in {1..22}; do
  CHR_INPUT=$PLINKDIR/chr${CHR}_files.txt
  CHR_OUTPUT=$PLINKDIR/chr${CHR}_files.merge.txt

  awk '{print $1".bed", $1".bim", $1".fam"}' "$CHR_INPUT" > "$CHR_OUTPUT"

  if [ "$CHR" -ne 1 ]; then
    cat "$CHR_OUTPUT" >> "$MERGE_LIST"
  fi
done

