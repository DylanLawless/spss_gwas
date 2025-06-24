#!/bin/bash
#SBATCH --job-name=plinkrelated
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=03:00:00
#SBATCH --mem=90gb
#SBATCH -o ./log/plinkrelated.%j.out
#SBATCH -e ./log/plinkrelated.%j.err

PLINK=/work/gr-fe/lawless/tool/plink_1.90/plink
PLINKDIR=/work/gr-fe/lawless/spss/gwas/plink
INPUT=$PLINKDIR/replaced.ID.bim/SPSS.QC.impute.allfiles
REMOVE_IDS=$PLINKDIR/relatedness/ID.exclude.relatedness.txt
OUTPUT=$PLINKDIR/replaced.ID.bim/SPSS.QC.impute.allfiles.kinged

$PLINK \
  --bfile "$INPUT" \
  --remove "$REMOVE_IDS" \
  --make-bed \
  --out "$OUTPUT" \
  --threads 32

