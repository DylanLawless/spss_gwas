#!/bin/bash

PLINKOUT=/mnt/data1/lawless/spss/gwas/shcs
COUNTDIR=$PLINKOUT/count_varID_length
OUTFILE=$COUNTDIR/All.varID_over79.sed.list

# Step 1: extract original variant IDs
for file in "$COUNTDIR"/*.varID_over79; do
  awk -F '=' '{gsub(/^[ \t]+|[ \t]+$/, "", $1); print $1}' "$file" > "$file.original"
done

# Step 2: generate shortened replacement IDs
for file in "$COUNTDIR"/*.varID_over79; do
  sed 's/:[^:]*/:ID/2g' "$file" > "$file.short"
done

# Step 3: build sed rules
> "$OUTFILE"
for file in "$COUNTDIR"/*.varID_over79; do
  paste -d' ' "$file.original" "$file.short" | awk '{print "s/" $1 "/" $2 "/g"}' >> "$OUTFILE"
done

