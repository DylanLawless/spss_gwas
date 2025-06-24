#!/bin/bash
#SBATCH --job-name=mergeCasesControls
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=02:00:00
#SBATCH --mem=30gb
#SBATCH -o ./log/9.merge.%j.out
#SBATCH -e ./log/9.merge.%j.err

PLINK=/home/lawless/tool/plink
OUTDIR=/mnt/data1/lawless/spss/merged
mkdir -p "$OUTDIR"

SPSS_BFILE=/mnt/data1/lawless/spss/plink/replaced.ID.bim/SPSS.QC.impute.allfiles.kinged.cleaned.filt
SHCS_BFILE=/mnt/data1/lawless/spss/gwas/shcs/replaced.ID.bim/SHCS.QC.impute.allfiles.kinged.cleaned.filt

# Attempt to merge directly
$PLINK \
  --bfile "$SPSS_BFILE" \
  --bmerge "${SHCS_BFILE}.bed" "${SHCS_BFILE}.bim" "${SHCS_BFILE}.fam" \
  --make-bed \
  --out "$OUTDIR/SPSS_SHCS.merged"

