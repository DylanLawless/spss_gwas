#!/bin/bash
#SBATCH --job-name=printid
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=00:30:00
#SBATCH --mem=60gb
#SBATCH -o ./log/4.5.printIDover79.%j.out
#SBATCH -e ./log/4.5.printIDover79.%j.err
#SBATCH --reservation=gr-fe

PLINKOUT=/mnt/data1/lawless/spss/gwas/shcs
COUNTDIR=$PLINKOUT/count_varID_length

for file in "$COUNTDIR"/SHCS.QC_chr*.count; do
  awk -F '=' '$2 > 79 { print $0 }' "$file" > "${file}.varID_over79"
done

