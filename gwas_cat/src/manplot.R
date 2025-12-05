#!/usr/bin/env Rscript

assoc_file <- "../data/cohort.filt.group.case_control.loco.mlma.printed.sig_snps_only_noSNP.tsv_formatted.tsv"
out_dir <- "../images/"

# Usage example from bash:
# Rscript manplot.R gwas_results.assoc.logistic /path/to/out_dir hard_coded

# args <- commandArgs(trailingOnly = TRUE)
# assoc_file <- args[1]
# out_dir <- args[2]
# run_label <- args[3]  # e.g., "hard_coded" or "qv_set_x"
run_label <- "spss"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ggmanh", quietly = TRUE)) BiocManager::install("ggmanh")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

suppressPackageStartupMessages({
  library(ggmanh)
  library(ggplot2)
})

# ---- read GWAS results ----
dat <- read.table(assoc_file, header = TRUE, stringsAsFactors = FALSE)
dat <- subset(dat, TEST == "ADD")
dat_hold <- dat
dat <- dat_hold
# dat <- head(dat)
dat$P <- suppressWarnings(as.numeric(dat$p_value))
dat <- dat[is.finite(dat$P) & dat$P > 0 & dat$P <= 1, ]

# ---- extract QC summary from log file ----
# log_file <- file.path(dirname(assoc_file), "gwas_results.log")
# qc_subtitle <- ""
# if (file.exists(log_file)) {
#   log_lines <- readLines(log_file, warn = FALSE)
#   keep <- grep(
#     "Before main variant filters|Total genotyping rate|pass filters and QC|cases and .*controls",
#     log_lines,
#     value = TRUE,
#     ignore.case = TRUE
#   )
#   clean <- gsub("^\\s*-*\\s*extract:\\s*", "", keep, ignore.case = TRUE)
#   clean <- gsub("^\\s*-+\\s*", "", clean)
#   clean <- gsub("\\s+", " ", trimws(clean))
#   wrapped <- vapply(clean, function(s) paste(strwrap(s, width = 80), collapse = "\n"), character(1))
#   qc_subtitle <- paste(wrapped, collapse = "\n")
# }

# ---- lambda ----
chisq <- qchisq(1 - dat$P, 1)
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
cat(sprintf("Genomic inflation factor (lambda): %.3f\n", lambda))

# ---- prepare df for ggmanh ----
# df <- data.frame(
#   chromosome = ifelse(dat$chromosome == 23, "X", as.character(dat$chromosome)),
#   position = dat$base_pair_location,
#   P = dat$p_value,
#   SNP = dat$variant_id,
#   stringsAsFactors = FALSE
# )

dat$SNP <- paste0(dat$chromosome, ":", dat$base_pair_location)

df <- data.frame(
  chromosome = ifelse(dat$chromosome == 23, "X", as.character(dat$chromosome)),
  position = dat$base_pair_location,
  P = dat$P,   # already present
  SNP = dat$SNP,
  stringsAsFactors = FALSE
)

df$chromosome <- factor(df$chromosome, c(as.character(1:22), "X"))

# ---- manhattan ----
g1 <- manhattan_plot(
  x = df,
  pval.colname = "P",
  chr.colname = "chromosome",
  pos.colname = "position",
  plot.title = sprintf("GWAS Manhattan (%s)", run_label),
  # plot.subtitle = qc_subtitle,
  y.label = "P",
  rescale = FALSE
)
ggsave(file.path(out_dir, sprintf("manhattan_%s.png", run_label)), plot = g1, width = 6, height = 4)

# ---- qq plot ----
qq_df <- within(df[is.finite(df$P), , drop = FALSE], {
  obs = -log10(sort(P))
  exp = -log10(ppoints(length(obs)))
})
qq_plot <- ggplot(qq_df, aes(x = exp, y = obs)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_point(size = 0.6, alpha = 0.7) +
  labs(
    title = sprintf("QQ plot  lambda = %.3f (%s)", lambda, run_label),
    # subtitle = qc_subtitle,
    x = "expected -log10(p)",
    y = "observed -log10(p)"
  ) +
  theme_bw()
ggsave(file.path(out_dir, sprintf("qqplot_%s.png", run_label)), qq_plot, width = 6, height = 6)

# ---- session info ----
# writeLines(c(capture.output(sessionInfo())), con = file.path(out_dir, sprintf("sessionInfo_%s.txt", run_label)))

