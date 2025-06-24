#!/bin/bash
#SBATCH --job-name=convertbfile
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=04:00:00
#SBATCH --mem=60gb
#SBATCH -o ./log/3.convertbfile.%j.out
#SBATCH -e ./log/3.convertbfile.%j.err
#SBATCH --reservation=gr-fe

PLINK=/work/gr-fe/lawless/tool/plink_1.90/plink
FILELIST_DIR=/mnt/data1/lawless/spss/gwas/filelist
PLINKOUT=/mnt/data1/lawless/spss/gwas/shcs

for CHR in {1..22}; do
  while read -r line; do
    $PLINK --file "$PLINKOUT/$line" --make-bed --out "$PLINKOUT/$line"
    $PLINK --file "$PLINKOUT/$line" --recode --out "$PLINKOUT/${line}.recode"
  done < "$FILELIST_DIR/chr${CHR}_files"
done

