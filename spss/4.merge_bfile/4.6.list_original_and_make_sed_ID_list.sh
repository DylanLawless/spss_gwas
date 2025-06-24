#!/bin/bash

PLINKOUT=/mnt/data1/lawless/spss/plink
COUNTDIR=$PLINKOUT/count_varID_length
OUTFILE=$COUNTDIR/All.varID_over79.sed.list

# [1] Extract original long IDs
for file in "$COUNTDIR"/*varID_over79; do
  awk -F '=' '{gsub(/^[ \t]+|[ \t]+$/, "", $1); print $1}' "$file" > "$file.original_longID"
done

# [2] Create shortened replacement IDs by modifying the second colon group
for file in "$COUNTDIR"/*varID_over79; do
  sed 's/:[^:]*/:ID/2g' "$file" > "$file.ID"
done

# [3] Paste original and modified into sed replacement format: s/old/new/g
for file in "$COUNTDIR"/*varID_over79; do
  paste -d' ' "$file.original_longID" "$file.ID" | \
    awk '{print "s/" $1 "/" $2 "/g"}' > "$file.sed.list"
done

# [4] Concatenate all into one sed script
cat "$COUNTDIR"/*.sed.list > "$OUTFILE"

