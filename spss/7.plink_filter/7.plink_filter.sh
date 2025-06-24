#!/bin/bash
#SBATCH --job-name=plinkfilt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=90gb
#SBATCH -o ./log/plinkfilt.%j.out
#SBATCH -e ./log/plinkfilt.%j.err

PLINK=/work/gr-fe/lawless/tool/plink_1.90/plink
PLINKDIR=/work/gr-fe/lawless/spss/gwas/plink
INPUT=$PLINKDIR/replaced.ID.bim/SPSS.QC.impute.allfiles.kinged.cleaned
OUTPUT=$PLINKDIR/replaced.ID.bim/SPSS.QC.impute.allfiles.kinged.cleaned.filt

$PLINK \
  --bfile "$INPUT" \
  --maf 0.01 \
  --geno 0.2 \
  --hwe 0.0000001 \
  --make-bed \
  --allow-no-sex \
  --out "$OUTPUT"

