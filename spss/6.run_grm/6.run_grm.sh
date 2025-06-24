#!/bin/bash
#SBATCH --job-name=grm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=03:00:00
#SBATCH --mem=90gb
#SBATCH -o ./log/grm.%j.out
#SBATCH -e ./log/grm.%j.err

PLINKDIR=/work/gr-fe/lawless/spss/gwas/plink
GRMOUTPUT=$PLINKDIR/grm
GCTA=/work/gr-fe/lawless/tool/gcta64
INPUT=$PLINKDIR/replaced.ID.bim/SPSS.QC.impute.allfiles.kinged

mkdir -p "$GRMOUTPUT"

# Step 1: Make GRM
$GCTA \
  --bfile "$INPUT" \
  --autosome \
  --make-grm \
  --out "$GRMOUTPUT/SPSS.QC.impute.allfiles.kinged.grm" \
  --thread-num 32

# Step 2: Run PCA
$GCTA \
  --grm "$GRMOUTPUT/SPSS.QC.impute.allfiles.kinged.grm" \
  --pca 20 \
  --out "$GRMOUTPUT/SPSS.QC.impute.allfiles.kinged.pca" \
  --thread-num 18

