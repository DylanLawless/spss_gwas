#!/bin/bash
#SBATCH --job-name=mlma_caseonly
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=12:00:00
#SBATCH --mem=90gb
#SBATCH -o ./log/mlma_caseonly.%j.out
#SBATCH -e ./log/mlma_caseonly.%j.err

GCTA=/work/gr-fe/lawless/tool/gcta64
PROJECT=/work/gr-fe/lawless/spss/gwas/plink
GENO=$PROJECT/replaced.ID.bim/SPSS.QC.impute.allfiles.kinged.cleaned.filt
GRM=$PROJECT/grm/SPSS.QC.impute.allfiles.kinged.grm
PHENO_DIR=$PROJECT/phenotypes
RESULTS=$PROJECT/results
COV_CAT=$PHENO_DIR/covariates_cat
COV_QUANT=$PHENO_DIR/covariates_quant
PHENO_LIST=$PHENO_DIR/pheno_list.txt

# pheno_list.txt
# pheno_comorbidity
# pheno_gn
# pheno_icu_binary
# pheno_psofa.score
# pheno_pelod.score
# pheno_hospital.acquired
# pheno_outcome.death
# pheno_cvc.clabsi
# pheno_risk.category

mkdir -p "$RESULTS"

while read PHENO; do
  $GCTA \
    --mlma-loco \
    --bfile "$GENO" \
    --grm "$GRM" \
    --pheno "$PHENO_DIR/$PHENO" \
    --covar "$COV_CAT" \
    --qcovar "$COV_QUANT" \
    --out "$RESULTS/SPSS.QC.impute.allfiles.kinged.cleaned.filt_${PHENO}" \
    --thread-num 28
done < "$PHENO_LIST"



