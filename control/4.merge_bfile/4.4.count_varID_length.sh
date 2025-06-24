#!/bin/bash
#SBATCH --job-name=countvarid
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=00:30:00
#SBATCH --mem=60gb
#SBATCH -o ./log/4.4.countVarID.%j.out
#SBATCH -e ./log/4.4.countVarID.%j.err
#SBATCH --reservation=gr-fe

FILELIST_DIR=/mnt/data1/lawless/spss/gwas/filelist
PLINKOUT=/mnt/data1/lawless/spss/gwas/shcs
OUTDIR=$PLINKOUT/count_varID_length

mkdir -p "$OUTDIR"

for CHR in {1..22}; do
  FILELIST="$FILELIST_DIR/chr${CHR}_files"
  while read -r line; do
    awk '{ print $2 " = " length($2) }' "$PLINKOUT/$line.bim" > "$OUTDIR/$line.count"
  done < "$FILELIST"
done

