library(qqman)
library(dplyr)

#==========================
# Manhattan plot (raw)
#==========================
results_log <- read.table("./cohort.filt.group.case_control.loco.mlma.printed.sig", header=TRUE) %>%
  filter(p > 0)

png("./images/5.case_control.png", width = 8, height = 5, units = 'in', res = 350)
manhattan(results_log, chr="Chr", bp="bp", p="p",
          suggestiveline = FALSE, genomewideline = FALSE)
dev.off()

#==========================
# QQ plot (raw)
#==========================
png("./images/6.case_control-QQ.png", width = 8, height = 8, units = 'in', res = 350)
qq(results_log$p)
dev.off()

#==========================
# Lambda
#==========================
pvalue <- read.table("./pval", header=TRUE, sep=" ", na.strings=c("", "NA"))
chisq <- qchisq(1 - as.numeric(pvalue$p), 1)
lambda <- median(chisq) / qchisq(0.5, 1)
cat(lambda, "\n")

#==========================
# Near-significant SNPs (for trimming)
#==========================
tmp <- results_log %>%
  filter((p < 5e-6 & Chr < 9) | (p < 5e-6 & Chr > 9)) %>%
  sample_n(150) %>%
  select(SNP)

write.table(tmp, file = "tmp", sep = " ", quote = FALSE, col.names = NA)

# Exclude sampled SNPs
results_log2 <- filter(results_log, !(SNP %in% tmp$SNP))

#==========================
# Manhattan plot (cleaned)
#==========================
results_log_clean <- read.table("./cohort.filt.group.case_control.loco.mlma.printed.sig.150trim", header=TRUE) %>%
  filter(p > 0)

png("./images/5.clean_case_control.png", width = 8, height = 5, units = 'in', res = 350)
manhattan(results_log_clean, chr="Chr", bp="bp", p="p",
          suggestiveline = FALSE, genomewideline = FALSE)
dev.off()

#==========================
# QQ plot (cleaned)
#==========================
png("./images/6.clean_case_control-QQ.png", width = 8, height = 8, units = 'in', res = 350)
qq(results_log_clean$p)
dev.off()

