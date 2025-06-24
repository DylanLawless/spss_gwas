#!/bin/sh
#SBATCH --job-name=KING
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=08:00:00
#SBATCH --mem=90gb
#SBATCH -o ./log/king.%j.out
#SBATCH -e ./log/king.%j.err

KING=/work/gr-fe/lawless/tool/king
INPUT=/work/gr-fe/lawless/spss/gwas/mnt/plink/replaced.ID.bim/SPSS.QC.impute.allfiles.bed

$KING -b "$INPUT" --kinship --degree 2 --cpus 32

