library(EnsDb.Hsapiens.v75)

ld_token <- "d780...."
# ld_token <- Sys.getenv("LDLINK_TOKEN")

library(tidyverse)
library(patchwork)

# read data
df <- readr::read_tsv(
  "print_significant_snp.txt",
  col_types = cols(
    Chr = col_character(),
    SNP = col_character(),
    bp = col_double(),
    A1 = col_character(),
    A2 = col_character(),
    Freq = col_double(),
    b = col_double(),
    se = col_double(),
    p = col_double()
  )
) %>%
  mutate(
    pos_mb = bp / 1e6,
    neglogp = -log10(p)
  )

# allele frequency plot
p_freq <- ggplot(df, aes(pos_mb, Freq)) +
  geom_point(size = 3, colour = "#1f77b4") +
  geom_line(aes(group = 1), linewidth = 0.4, colour = "#1f77b4", alpha = 0.6) +
  geom_text(aes(label = SNP), nudge_y = 0.03, size = 3, check_overlap = TRUE) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "position on chr9 (Mb)", y = "allele frequency (A1)", title = "Allele frequency") +
  theme_bw(base_size = 12)

# effect size with SE
p_beta <- ggplot(df, aes(pos_mb, b)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_pointrange(aes(ymin = b - se, ymax = b + se), colour = "#2ca02c", fatten = 1, linewidth = 0.4) +
  labs(x = "position on chr9 (Mb)", y = "effect size (beta)", title = "Effect size with SE") +
  theme_bw(base_size = 12)

# standard error across position
p_se <- ggplot(df, aes(pos_mb, se)) +
  geom_point(size = 3, colour = "#9467bd") +
  geom_line(aes(group = 1), linewidth = 0.4, colour = "#9467bd", alpha = 0.6) +
  labs(x = "position on chr9 (Mb)", y = "standard error", title = "SE") +
  theme_bw(base_size = 12)

# -log10(p) across position
p_p <- ggplot(df, aes(pos_mb, neglogp)) +
  geom_point(size = 3, colour = "#d62728") +
  geom_line(aes(group = 1), linewidth = 0.4, colour = "#d62728", alpha = 0.6) +
  geom_hline(yintercept = 7.3, linetype = "dashed", colour = "grey60") +
  labs(x = "position on chr9 (Mb)", y = expression(-log[10](italic(p))), title = "-log10(p)") +
  theme_bw(base_size = 12)

# arrange
(p_p / p_freq / p_beta / p_se )



# Overlay  ----
library(tidyverse)
library(scales)
library(ggnewscale)

# read data
df <- readr::read_tsv(
  "print_significant_snp.txt",
  col_types = cols(
    Chr = col_character(),
    SNP = col_character(),
    bp = col_double(),
    A1 = col_character(),
    A2 = col_character(),
    Freq = col_double(),
    b = col_double(),
    se = col_double(),
    p = col_double()
  )
) %>%
  mutate(
    pos_mb = bp / 1e6,
    neglogp = -log10(p)
  )

# ranges for scaling
f_min <- min(df$Freq, na.rm = TRUE)
f_max <- max(df$Freq, na.rm = TRUE)
p_min <- min(df$neglogp, na.rm = TRUE)
p_max <- max(df$neglogp, na.rm = TRUE)

# scale -log10(p) to frequency scale for plotting
df <- df %>%
  mutate(p_scaled = rescale(neglogp, to = c(f_min, f_max), from = c(p_min, p_max)))

# plot with p-coloured connectors and separate colour scale for points
ggplot(df, aes(x = pos_mb)) +
  geom_segment(
    aes(x = pos_mb, xend = pos_mb, y = Freq, yend = p_scaled, colour = neglogp),
    linewidth = 0.7, alpha = 0.6
  ) +
  scale_colour_viridis_c(
    option = "plasma",
    name = expression(-log[10](italic(p)))
  ) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(y = Freq, colour = "allele frequency"), size = 3, alpha = 0.6) +
  geom_point(aes(y = p_scaled, colour = "-log10(p)"), size = 3, alpha = 0.6) +
  scale_colour_manual(
    name = NULL,
    values = c("allele frequency" = "#1f77b4", "-log10(p)" = "#d62728")
  ) +
  scale_y_continuous(
    name = "allele frequency (A1)",
    # limits = c(0, 1),
    sec.axis = sec_axis(
      trans = ~ rescale(., to = c(p_min, p_max), from = c(f_min, f_max)),
      name = expression(-log[10](italic(p)))
    )
  ) +
  labs(x = "position on chr9 (Mb)", title = "allele frequency and -log10(p) across region") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")


# p1 locus ----

# minimal locuszoom-style plot using gene + flank, with LD colouring if available

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(scales)
  library(ggplot2); ggplot2::theme_set(ggplot2::theme_bw())
  library(locuszoomr)
})

# example data in memory:
# data(SLE_gwas_sub)

# or read your own GWAS table with columns: chrom, pos, rsid, p [, r2]
gwas <- readr::read_tsv("print_significant_snp.txt", show_col_types = FALSE) |>
  dplyr::transmute(chrom = Chr, pos = bp, rsid = SNP, p = p)

loc <- locus(data = gwas, gene = 'FAM206A', flank = 1e5,
             ens_db = "EnsDb.Hsapiens.v75")
summary(loc)
locus_plot(loc)



# gwas <- SLE_gwas_sub |>
  # dplyr::mutate(chrom = as.character(chrom)) |>
  # dplyr::select(chrom, pos, rsid, p, dplyr::any_of("r2"))

gene_symbol <- "FAM206A"
flank_bp <- 1e5
ensdb_pkg <- "EnsDb.Hsapiens.v75"  # use "EnsDb.Hsapiens.v86" for GRCh38

# 1) build locus object using gene + flank
loc <- locuszoomr::locus(
  data   = gwas,
  gene   = gene_symbol,
  flank  = flank_bp,
  ens_db = ensdb_pkg
)

# 2) map colours
if ("r2" %in% names(gwas) && any(is.finite(gwas$r2))) {
  ld_bins <- c("<0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", ">=0.8")
  ld_cols <- c("<0.2" = "#4575b4", "0.2–0.4" = "#91bfdb", "0.4–0.6" = "#fee090",
               "0.6–0.8" = "#fc8d59", ">=0.8" = "#d73027")
  gwas <- gwas |>
    dplyr::mutate(
      ld_bin = cut(r2, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
                   labels = ld_bins, right = FALSE),
      bg = unname(ld_cols[as.character(ld_bin)])
    )
} else {
  gwas <- gwas |>
    dplyr::mutate(
      neglogp = -log10(p),
      bg = viridisLite::plasma(100)[scales::rescale(neglogp, to = c(1, 100)) |> round()]
    )
}
gwas$bg <- scales::alpha(gwas$bg, 0.85)

# 3) join colours into loc$data
loc$data <- loc$data |>
  dplyr::left_join(gwas |> dplyr::select(rsid, bg), by = "rsid")

bg_vec <- ifelse(is.na(loc$data$bg), "#99999980", loc$data$bg)

# 4) plot
locuszoomr::locus_plot(
  loc,
  bg = bg_vec,
  ylab = expression(-log[10](italic(p))),
  pch = 21,
  col = "black",
  cex = 1.2,
  pcutoff = NULL,
  cex.axis = 0.85,
  cex.lab = 1.0,
  cex.main = 1.0
)



# p2 locus ----

# locuszoom-style plot using LD from LDlinkR

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(scales)
  library(ggplot2); ggplot2::theme_set(ggplot2::theme_bw())
  library(locuszoomr)
})

# config
gene_symbol <- "FAM206A"
flank_bp <- 1e5
ensdb_pkg <- "EnsDb.Hsapiens.v75"  # use "EnsDb.Hsapiens.v86" for GRCh38
ld_pop <- "CEU" #  (e.g. EUR, YRI or CEU),
ld_build <- "GRCh37"               # set "GRCh38" if needed
ld_build <- "grch37"
# ld_token <- Sys.getenv("LDLINK_TOKEN")

# read GWAS table with columns: Chr, bp, SNP, p
gwas <- readr::read_tsv("print_significant_snp.txt", show_col_types = FALSE) |>
  dplyr::transmute(
    chrom = as.character(Chr),
    pos = as.numeric(bp),
    rsid = as.character(SNP),
    p = as.numeric(p),
    neglogp = -log10(p)
  ) |>
  dplyr::filter(is.finite(neglogp), !is.na(rsid), !is.na(pos))

# build locus object using gene and flank to define window
loc <- locuszoomr::locus(
  data   = gwas |> dplyr::select(chrom, pos, rsid, p),
  gene   = gene_symbol,
  flank  = flank_bp,
  ens_db = ensdb_pkg
)

# filter gwas to the window used by locus
loc_chr <- as.character(loc$data$chrom)
loc_start <- as.numeric(loc$TX$start)
loc_end <- as.numeric(loc$TX$end)

region <- gwas |>
  dplyr::filter(
    chrom == loc_chr,
    pos >= loc_start,
    pos <= loc_end
  ) |>
  dplyr::arrange(pos)

# choose lead SNP as the smallest p in the window
lead_rsid <- region |> dplyr::arrange(p) |> dplyr::slice(1) |> dplyr::pull(rsid)
lead_rsid <- region |> dplyr::arrange(p) |> dplyr::slice(17) |> dplyr::pull(rsid)


# fetch LD matrix from LDlinkR
if (!requireNamespace("LDlinkR", quietly = TRUE)) {
  stop("Package LDlinkR is required to download LD. Please install LDlinkR.")
}
if (!nzchar(ld_token)) {
  stop("LDLINK_TOKEN is missing. Set Sys.setenv(LDLINK_TOKEN = 'your_token') before running.")
}

rs_list <- unique(region$rsid)
# LDlink recommends fewer than about 300 SNPs per call. If more, keep the nearest 300 by position
if (length(rs_list) > 300) {
  idx <- order(abs(region$pos - region$pos[match(lead_rsid, region$rsid)]))
  rs_list <- unique(region$rsid[idx])[seq_len(300)]
}

ldmat <- LDlinkR::LDmatrix(
  snps = rs_list,
  pop = ld_pop,
  r2d = "r2",
  genome_build = ld_build,
  token = ld_token,
  file = FALSE
)

# extract r2 to lead and join back
rs_col <- dplyr::coalesce(ldmat$RS_number, ldmat[[1]])
rownames(ldmat) <- rs_col
colnames(ldmat)[1] <- "RS_number"

if (!(lead_rsid %in% colnames(ldmat))) {
  stop(paste0("Lead SNP ", lead_rsid, " not present in LD matrix."))
}

r2_vec <- suppressWarnings(as.numeric(unlist(ldmat[, lead_rsid])))
names(r2_vec) <- rs_col

region <- region |>
  dplyr::left_join(
    tibble::tibble(rsid = names(r2_vec), r2 = as.numeric(r2_vec)),
    by = "rsid"
  ) |>
  dplyr::mutate(r2 = dplyr::if_else(rsid == lead_rsid & is.na(r2), 1, r2))

# colour strictly by LD bins
ld_bins <- c("<0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", ">=0.8")
ld_cols <- c(
  "<0.2" = "#4575b4",
  "0.2–0.4" = "#91bfdb",
  "0.4–0.6" = "#fee090",
  "0.6–0.8" = "#fc8d59",
  ">=0.8" = "#d73027"
)

region <- region |>
  dplyr::mutate(
    ld_bin = cut(r2, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
                 labels = ld_bins, right = FALSE),
    bg = scales::alpha(unname(ld_cols[as.character(ld_bin)]), 0.85)
  )

# join colours into loc$data
loc$data <- loc$data |>
  dplyr::left_join(region |> dplyr::select(rsid, bg, r2, ld_bin), by = "rsid")

bg_vec <- ifelse(is.na(loc$data$bg), "#99999980", loc$data$bg)

# association plot coloured by LD
locuszoomr::locus_plot(
  loc,
  bg = bg_vec,
  ylab = expression(-log[10](italic(p))),
  pch = 21,
  col = "black",
  cex = 1.2,
  pcutoff = NULL,
  cex.axis = 0.9,
  cex.lab = 1.0,
  cex.main = 1.0
)

# LD on Y confirmation plot
ggplot2::ggplot(loc$data, ggplot2::aes(x = pos / 1e6, y = r2)) +
  ggplot2::geom_point(
    ggplot2::aes(fill = ld_bin),
    colour = "black",
    shape = 21,
    size = 2.2,
    alpha = 0.85
  ) +
  ggplot2::scale_fill_manual(
    values = ld_cols,
    na.value = "#99999980",
    name = expression(LD~r^2~to~lead)
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::geom_hline(yintercept = c(0.2, 0.4, 0.6, 0.8), colour = "grey85", linewidth = 0.25) +
  ggplot2::labs(
    x = paste0("position on chr", loc_chr <- loc$region$chrom, " (Mb)"),
    y = expression(LD~r^2),
    title = paste0("LD to lead variant ", lead_rsid, " across ", gene_symbol)
  ) +
  ggplot2::theme(legend.position = "right")





# 1) fetch proxies around the lead SNP to expand the set
library(dplyr)
library(LDlinkR)

lead <- lead_rsid
prox <- LDproxy(
  snp = lead,
  pop = ld_pop,              # "EUR", etc.
  r2d = "r2",
  token = ld_token,
  genome_build = "grch37"   # exact spelling
  # output = "table"
)

# prox has cols RS_Number, Coordinate, R2, etc.
prox_tbl <- prox |>
  transmute(
    rsid = RS_Number,
    r2 = as.numeric(R2)
  ) |>
  distinct(rsid, .keep_all = TRUE)

# 2) merge proxies with your GWAS table
region_expanded <- region |>
  full_join(prox_tbl, by = "rsid")

# keep p from GWAS when present, set missing p to 1 so they sit at y=0
region_expanded <- region_expanded |>
  mutate(
    p = if_else(is.na(p), 1, p),
    neglogp = -log10(p)
  )

# colour strictly by LD and replot

# 3) LD bins and colours
ld_cols <- c(
  "<0.2" = "#4575b4",
  "0.2–0.4" = "#91bfdb",
  "0.4–0.6" = "#fee090",
  "0.6–0.8" = "#fc8d59",
  ">=0.8" = "#d73027"
)

region_expanded <- region_expanded |>
  dplyr::mutate(
    ld_bin = cut(
      r2.x,
      breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
      labels = names(ld_cols),
      right = FALSE
    ),
    bg = scales::alpha(unname(ld_cols[as.character(ld_bin)]), 0.85)
  )

# 4) association plot coloured by LD
p_assoc_ld <- ggplot2::ggplot(
  region_expanded,
  ggplot2::aes(x = pos / 1e6, y = neglogp)
) +
  ggplot2::geom_point(
    ggplot2::aes(fill = ld_bin),
    colour = "black",
    shape = 21,
    size = 2.2,
    alpha = 0.85
  ) +
  ggplot2::scale_fill_manual(values = ld_cols, na.value = "#99999980", name = expression(LD~r^2~to~lead)) +
  ggplot2::labs(x = "position (Mb)", y = expression(-log[10](italic(p)))) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "top")

print(p_assoc_ld)

# 5) LD on Y confirmation plot
p_ld_y <- ggplot2::ggplot(
  region_expanded,
  ggplot2::aes(x = pos / 1e6, y = r2.x)
) +
  ggplot2::geom_point(
    ggplot2::aes(fill = ld_bin),
    colour = "black",
    shape = 21,
    size = 2.2,
    alpha = 0.85
  ) +
  ggplot2::scale_fill_manual(values = ld_cols, na.value = "#99999980", name = expression(LD~r^2~to~lead)) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::geom_hline(yintercept = c(0.2, 0.4, 0.6, 0.8), colour = "grey85", linewidth = 0.25) +
  ggplot2::labs(x = "position (Mb)", y = expression(LD~r^2)) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "top")

print(p_ld_y)

















# independence -----
# approach 1: fine-map with susie_rss using public LD to infer independent signals

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(LDlinkR)
  library(susieR)
  library(Matrix)
})

# inputs you already have
ld_pop   <- "EUR"        # choose population close to Swiss cohort
ld_build <- "grch37"     # exact spelling for LDlinkR
ld_token <- ld_token

gwas <- read_tsv("print_significant_snp.txt", show_col_types = FALSE) |>
  transmute(
    rsid = as.character(SNP),
    pos  = as.numeric(bp),
    beta = as.numeric(b),
    se   = as.numeric(se),
    p    = as.numeric(p)
  ) #|>
  # filter(is.finite(beta), is.finite(se), !is.na(rsid))

# choose region and lead
lead_rsid <- gwas |> arrange(p) |> dplyr::slice(1) |> pull(rsid)
rs_list    <- unique(gwas$rsid)

stopifnot(requireNamespace("LDlinkR", quietly = TRUE), nzchar(ld_token))

# download LD matrix among your rsids
ldmat <- LDmatrix(
  snps = rs_list,
  pop = ld_pop,
  r2d = "r2",
  token = ld_token,
  genome_build = ld_build,
  file = FALSE
)

# tidy LD matrix to numeric symmetric R
rs_col <- dplyr::coalesce(ldmat$RS_number, ldmat[[1]])
rownames(ldmat) <- rs_col
colnames(ldmat)[1] <- "RS_number"
R <- as.matrix(ldmat[,-1, drop = FALSE])
rownames(R) <- rs_col
colnames(R) <- colnames(R)

# harmonise order with GWAS z-scores
g <- gwas |> dplyr::filter(rsid %in% rownames(R))
R <- R[g$rsid, g$rsid, drop = FALSE]
R <- suppressWarnings(nearPD(R, corr = TRUE)$mat) |> as.matrix()

# z-scores and sample size
z <- g$beta / g$se
# provide your study N if known, else a proxy
N <- 510 + 994

# run SuSiE in LD space
fit <- susie_rss(z = z, R = R, n = N, L = 3, max_iter = 200)

# extract credible sets and lead variants per set
cs <- susie_get_cs(fit, X = NULL)
signals <- tibble(
  set = seq_along(cs$cs),
  cs_size = lengths(cs$cs),
  lead_index = sapply(cs$cs, function(idx) idx[which.max(fit$pip[idx])]),
  lead_rsid = g$rsid[unlist(lapply(cs$cs, function(idx) idx[which.max(fit$pip[idx])]))],
  lead_pip  = sapply(cs$cs, function(idx) max(fit$pip[idx]))
)

signals






# independence -----
# approach 1: fine-map with susie_rss using public LD to infer independent signals
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(LDlinkR)
  library(susieR)
  library(Matrix)
})

cat("Starting fine-mapping with public LD\n")

# inputs you already have
ld_pop   <- "EUR"         # population closest to Swiss cohort
ld_build <- "grch37"      # genome build for LDlinkR, as provided
ld_token <- ld_token      # assumes your token object is already set

cat(paste0("LD population: ", ld_pop, "\n"))
cat(paste0("LD genome build: ", ld_build, "\n"))

# read GWAS
gwas <- readr::read_tsv("print_significant_snp.txt", show_col_types = FALSE) |>
  dplyr::transmute(
    rsid = as.character(.data$SNP),
    pos  = as.numeric(.data$bp),
    beta = as.numeric(.data$b),
    se   = as.numeric(.data$se),
    p    = as.numeric(.data$p)
  )

cat(paste0("Total variants read: ", nrow(gwas), "\n"))

# choose region and lead
lead_rsid <- gwas |> dplyr::arrange(.data$p) |> dplyr::slice(1) |> dplyr::pull(.data$rsid)
rs_list   <- unique(gwas$rsid)

cat(paste0("Lead variant by P: ", lead_rsid, "\n"))
cat(paste0("Requesting LD for ", length(rs_list), " rsIDs\n"))

stopifnot(requireNamespace("LDlinkR", quietly = TRUE), nzchar(ld_token))

# download LD matrix among your rsids
# ldmat <- LDlinkR::LDmatrix(
#   snps = rs_list,
#   pop = ld_pop,
#   r2d = "r2",
#   token = ld_token,
#   genome_build = ld_build,
#   file = FALSE
# )



# cache LD matrix to disk and only download if missing

# set up a simple cache path
cache_dir <- "ld_cache"
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# build a cache key from population, build, and the SNP list size
# keep filename short and filesystem-safe
cache_key <- paste0(
  "LDmatrix_", ld_pop, "_", ld_build, "_n", length(rs_list), ".rds"
)
cache_path <- file.path(cache_dir, cache_key)

if (file.exists(cache_path)) {
  cat("Reading LD matrix from cache: ", cache_path, "\n", sep = "")
  ldmat <- readRDS(cache_path)
} else {
  cat("Downloading LD matrix from LDlinkR\n")
  ldmat <- LDlinkR::LDmatrix(
    snps = rs_list,
    pop = ld_pop,
    r2d = "r2",
    token = ld_token,
    genome_build = ld_build,
    file = FALSE
  )
  saveRDS(ldmat, cache_path)
  cat("Saved LD matrix to: ", cache_path, "\n", sep = "")
}





# tidy LD matrix to numeric symmetric R
rs_col <- dplyr::coalesce(ldmat$RS_number, ldmat[[1]])
rownames(ldmat) <- rs_col
colnames(ldmat)[1] <- "RS_number"
R <- as.matrix(ldmat[ , -1, drop = FALSE])
rownames(R) <- rs_col
colnames(R) <- colnames(R)

cat(paste0("LD matrix received with ", nrow(R), " SNPs\n"))

# harmonise order with GWAS z-scores
g <- gwas |> dplyr::filter(.data$rsid %in% rownames(R))
R <- R[g$rsid, g$rsid, drop = FALSE]

cat(paste0("Variants overlapping between GWAS and LD: ", nrow(g), "\n"))

# ensure correlation matrix properties
R <- suppressWarnings(Matrix::nearPD(R, corr = TRUE)$mat) |> as.matrix()
diag_dev <- max(abs(diag(R) - 1))
symm_dev <- max(abs(R - t(R)))
cat(paste0("Max deviation of diag(R) from 1: ", signif(diag_dev, 3), "\n"))
cat(paste0("Max asymmetry of R: ", signif(symm_dev, 3), "\n"))

# z-scores and sample size
z <- g$beta / g$se
N <- 510 + 994

cat(paste0("Computed z-scores for ", length(z), " variants\n"))
cat(paste0("Using sample size N = ", N, "\n"))

# run SuSiE in LD space
cat("Running SuSiE fine-mapping\n")
fit <- susieR::susie_rss(z = z, R = R, n = N, L = 3, max_iter = 200, coverage = 0.95, verbose = FALSE)

cat(paste0("SuSiE converged: ", ifelse(isTRUE(fit$converged), "yes", "no"), "\n"))

# extract credible sets and lead variants per set
cs <- susieR::susie_get_cs(fit, X = NULL)

n_sets <- length(cs$cs)
cat(paste0("Number of 95% credible sets: ", n_sets, "\n"))

if (n_sets == 0) {
  cat("Interpretation: no distinct credible set was detected given the public LD. The region likely reflects a single weak signal or insufficient information.\n")
} else {
  pip <- fit$pip
  for (i in seq_len(n_sets)) {
    idx <- cs$cs[[i]]
    lead_index <- idx[which.max(pip[idx])]
    lead_rsid_cs <- g$rsid[lead_index]
    lead_pip  <- max(pip[idx])
    cs_size   <- length(idx)
    cs_cov    <- if (!is.null(cs$coverage)) cs$coverage[i] else NA_real_
    best_p    <- min(g$p[idx], na.rm = TRUE)
    
    cat(paste0(
      "Set ", i, ": size ", cs_size,
      ", coverage ", ifelse(is.na(cs_cov), "NA", sprintf("%.2f", cs_cov)),
      ", lead ", lead_rsid_cs,
      ", lead PIP ", sprintf("%.3f", lead_pip),
      ", best P ", format(best_p, scientific = TRUE), "\n"
    ))
    
    # print top variants in the set
    ord <- order(pip[idx], decreasing = TRUE)
    show_k <- min(5, length(ord))
    top <- data.frame(
      rank = seq_len(show_k),
      rsid = g$rsid[idx[ord[seq_len(show_k)]]],
      pip  = signif(pip[idx[ord[seq_len(show_k)]]], 3),
      p    = signif(g$p[idx[ord[seq_len(show_k)]]], 3),
      pos  = g$pos[idx[ord[seq_len(show_k)]]]
    )
    print(top, row.names = FALSE)
  }
  
  if (n_sets == 1) {
    cat("Interpretation: one credible set suggests a single independent association under the public LD used.\n")
  } else {
    cat("Interpretation: multiple credible sets suggest more than one independent association under the public LD used.\n")
  }
}

cat("Notes:\n")
cat("- This analysis uses public LD and your GWAS summary effects. Swiss-specific LD may differ.\n")
cat("- Credible sets reflect statistical independence after conditioning on the LD structure and z-scores.\n")
cat("- If sets are large or unstable, refine variant inclusion or try a population panel closest to your cohort.\n")





# plot SuSiE-labelled variants with annotations ----

suppressPackageStartupMessages({
  library(ggplot2); ggplot2::theme_set(ggplot2::theme_bw())
  library(ggrepel)
  library(scales)
  library(dplyr)
})

# build a data frame with SuSiE results aligned to GWAS rows
susie_df <- g %>%
  dplyr::mutate(
    neglogp = -log10(p),
    pip = fit$pip[match(rsid, g$rsid)],
    cs_id = NA_integer_
  )

# assign credible set ids
if (length(cs$cs) > 0) {
  for (i in seq_along(cs$cs)) {
    idx <- cs$cs[[i]]
    susie_df$cs_id[idx] <- i
  }
}

# flag lead variant per credible set
lead_tbl <- NULL
if (length(cs$cs) > 0) {
  lead_tbl <- lapply(seq_along(cs$cs), function(i) {
    idx <- cs$cs[[i]]
    pip_i <- fit$pip[idx]
    lead_idx <- idx[which.max(pip_i)]
    data.frame(
      rsid = g$rsid[lead_idx],
      cs_id = i,
      lead = TRUE,
      stringsAsFactors = FALSE
    )
  }) %>% dplyr::bind_rows()
}

susie_df <- susie_df %>%
  dplyr::left_join(lead_tbl, by = c("rsid", "cs_id")) %>%
  dplyr::mutate(
    lead = dplyr::if_else(is.na(lead), FALSE, lead),
    cs_fac = factor(cs_id, levels = sort(unique(cs_id)))
  )

# label text for lead variants
label_df <- susie_df %>%
  dplyr::filter(lead) %>%
  dplyr::mutate(
    label = paste0(
      "set ", cs_id,
      "\n", rsid,
      "\nPIP=", scales::number(pip, accuracy = 0.001),
      "  P=", format(p, digits = 3, scientific = TRUE)
    )
  )

# palette for credible sets
n_sets <- length(unique(na.omit(susie_df$cs_id)))
pal <- scales::hue_pal()(max(1, n_sets))
names(pal) <- levels(susie_df$cs_fac)

# main locus-style plot: position vs -log10(P), coloured by credible set, sized by PIP
p_locus <- ggplot2::ggplot(susie_df, ggplot2::aes(x = pos / 1e6, y = neglogp)) +
  ggplot2::geom_point(
    ggplot2::aes(colour = cs_fac, size = pip),
    alpha = 0.85,
    stroke = 0.2
  ) +
  ggplot2::scale_colour_manual(
    values = pal,
    na.value = "#99999966",
    name = "credible set"
  ) +
  ggplot2::scale_size(range = c(1.2, 4), name = "PIP") +
  ggrepel::geom_text_repel(
    data = label_df,
    ggplot2::aes(label = label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.size = 0.25,
    colour = "black",
    min.segment.length = 0
  ) +
  ggplot2::labs(
    x = "position (Mb)",
    y = expression(-log[10](italic(p))),
    title = "regional association with SuSiE credible sets and PIP"
  ) +
  ggplot2::theme(
    legend.position = "right"
  )

print(p_locus)

ggsave("p_locus.pdf", plot = p_locus, height = 3, width = 6)

# optional diagnostic: PIP vs -log10(P) to see alignment between evidence types
p_diag <- ggplot2::ggplot(susie_df, ggplot2::aes(x = neglogp, y = pip)) +
  ggplot2::geom_point(ggplot2::aes(colour = cs_fac), alpha = 0.85, size = 2) +
  ggplot2::scale_colour_manual(values = pal, na.value = "#99999966", name = "credible set") +
  ggplot2::labs(x = expression(-log[10](italic(p))), y = "PIP", title = "diagnostic: PIP versus -log10(P)") +
  ggplot2::theme(legend.position = "right")

print(p_diag)



# view ----


# LD-based haplotype blocks and visualisations from SuSiE inputs

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2); ggplot2::theme_set(ggplot2::theme_bw())
  library(ggrepel)
  library(scales)
  library(igraph)
})

# assume g (GWAS table used for SuSiE) and R (LD matrix used in SuSiE) exist

# 1) get correlation matrix R_corr and r2
R_corr <- if (max(R, na.rm = TRUE) > 1) sqrt(R) else R
R_corr[is.na(R_corr)] <- 0
R_corr <- (R_corr + t(R_corr)) / 2
diag(R_corr) <- 1
R2 <- R_corr^2

# ensure order by position
ord <- order(g$pos)
g_ord <- g[ord, , drop = FALSE]
R2_ord <- R2[ord, ord, drop = FALSE]

# 2) compute LD blocks via graph components at threshold
ld_thr <- 0.6
edges <- which(R2_ord >= ld_thr & upper.tri(R2_ord), arr.ind = TRUE)
g_graph <- if (nrow(edges) > 0) {
  igraph::graph_from_edgelist(
    cbind(g_ord$rsid[edges[, 1]], g_ord$rsid[edges[, 2]]),
    directed = FALSE
  )
} else {
  igraph::make_empty_graph(n = 0)
}

block_ids <- rep(NA_integer_, nrow(g_ord))
if (igraph::gorder(g_graph) > 0) {
  comp <- igraph::components(g_graph)$membership
  block_ids[match(names(comp), g_ord$rsid)] <- comp
}
g_ord$ld_block <- block_ids

# 3) r2 to each SuSiE lead per credible set
susie_leads <- NULL
if (exists("cs") && length(cs$cs) > 0) {
  susie_leads <- lapply(seq_along(cs$cs), function(i) {
    idx <- cs$cs[[i]]
    pip_i <- fit$pip[idx]
    lead_idx <- idx[which.max(pip_i)]
    data.frame(cs_id = i, rsid = g$rsid[lead_idx], stringsAsFactors = FALSE)
  }) |> dplyr::bind_rows()
}

# map r2 to each lead
r2_lead_long <- NULL
if (!is.null(susie_leads)) {
  r2_lead_long <- susie_leads |>
    dplyr::rowwise() |>
    dplyr::mutate(
      r2_vec = list(R2[g$rsid, rsid])
    ) |>
    dplyr::ungroup() |>
    tidyr::unnest_wider(r2_vec, names_sep = "_") |>
    tidyr::pivot_longer(
      cols = starts_with("r2_vec_"),
      names_to = "target",
      values_to = "r2_to_lead"
    ) |>
    dplyr::mutate(target = sub("^r2_vec_", "", target)) |>
    dplyr::rename(lead_rsid = rsid)
  # join r2_to_lead back to ordered GWAS
  g_ord <- g_ord |>
    dplyr::left_join(
      r2_lead_long |>
        dplyr::rename(rsid = target),
      by = "rsid"
    )
}

# 4) plot A: regional association coloured by LD block
pal_blocks <- scales::hue_pal()(max(1, length(unique(na.omit(g_ord$ld_block)))))
names(pal_blocks) <- sort(unique(na.omit(g_ord$ld_block)))

p_blocks <- ggplot2::ggplot(g_ord, ggplot2::aes(x = pos / 1e6, y = -log10(p))) +
  ggplot2::geom_point(
    ggplot2::aes(colour = factor(ld_block)),
    size = 2,
    alpha = 0.9
  ) +
  ggplot2::scale_colour_manual(values = pal_blocks, na.value = "#99999966", name = "LD block\n(r2 ≥ 0.6)") +
  ggplot2::labs(x = "position (Mb)", y = expression(-log[10](italic(p))), title = "regional association coloured by LD blocks") +
  ggplot2::theme(legend.position = "right")

print(p_blocks)



# 5) plot B: r2 to each SuSiE lead, faceted by credible set

# build a clean facet label in the main table
g_ord_fac <- g_ord %>%
  dplyr::mutate(
    facet = dplyr::if_else(!is.na(cs_id),
                           paste0("set ", cs_id, " lead ", lead_rsid),
                           NA_character_)
  )

# label rows are the lead variants per credible set
lead_labels <- g_ord_fac %>%
  dplyr::filter(!is.na(cs_id), rsid == lead_rsid) %>%
  dplyr::mutate(
    label = paste0(
      "set ", cs_id, "\n", rsid,
      "\nPIP=", scales::number(fit$pip[match(rsid, g$rsid)], accuracy = 0.001),
      "  P=", format(p, digits = 3, scientific = TRUE)
    )
  )

p_r2_facets <- ggplot2::ggplot(
  g_ord_fac %>% dplyr::filter(!is.na(r2_to_lead), !is.na(facet)),
  ggplot2::aes(x = pos / 1e6, y = -log10(p))
) +
  ggplot2::geom_point(
    ggplot2::aes(colour = r2_to_lead),
    size = 2.2
  ) +
  ggplot2::scale_colour_viridis_c(name = expression(r^2~to~lead), limits = c(0, 1)) +
  ggrepel::geom_text_repel(
    data = lead_labels,
    ggplot2::aes(x = pos / 1e6, y = -log10(p), label = label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.25
  ) +
  ggplot2::facet_wrap(~ facet, scales = "free_x") +
  ggplot2::labs(
    x = "position (Mb)",
    y = expression(-log[10](italic(p))),
    title = "association coloured by r2 to SuSiE lead"
  ) +
  ggplot2::theme(legend.position = "right")

print(p_r2_facets)

# compare Swiss z to reference LD expectations and show conditional residuals

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2); theme_set(theme_bw())
  library(ggrepel)
})

# g and R are the GWAS table and reference LD matrix used for SuSiE
# assume g has columns: rsid, pos, beta, se, p
# R is a correlation matrix ordered as g$rsid

z <- g$beta / g$se
names(z) <- g$rsid

# pick lead as smallest P in g
lead <- g$rsid[which.min(g$p)]
z_lead <- z[lead]
r_to_lead <- R[, lead]              # reference correlations to lead
r2_to_lead <- r_to_lead^2

df <- g %>%
  mutate(
    neglogp = -log10(p),
    r = r_to_lead[rsid],
    r2 = r^2,
    z_pred = r * z_lead,
    z_ratio = z / z_lead           # should track r if LD matches
  )

# 1) z_ratio vs r: if reference LD matches cohort, points sit on the 1:1 line
p_conditional_1 <- ggplot(df, aes(x = r, y = z_ratio)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_point(aes(colour = r2), alpha = 0.8) +
  scale_colour_viridis_c(name = expression(r^2~ref)) +
  labs(x = "reference r to lead", y = "z_i / z_lead", title = "cohort z versus reference LD") +
  theme(legend.position = "right")

# 2) conditional residuals using reference LD: z_cond = z - r * z_lead
# if one true signal, residuals should be near 0 across the window
df$z_cond <- z - df$r * z_lead
df$neglogp_cond <- -log10(2*pnorm(-abs(df$z_cond)))  # two-sided P from residual z

p_conditional_2 <- ggplot(df, aes(x = pos/1e6)) +
  geom_point(aes(y = neglogp, colour = "raw"), alpha = 0.85) +
  geom_point(aes(y = neglogp_cond, colour = "conditional\non lead (ref LD)"), alpha = 0.85) +
  scale_colour_manual(values = c("raw" = "#1f77b4", "conditional\non lead (ref LD)" = "#d62728"), name = NULL) +
  labs(x = "position (Mb)", y = expression(-log[10](italic(p))),
       title = "raw versus conditional association using reference LD") +
  theme(legend.position = "right")

print(p_conditional_1)
print(p_conditional_2)


ggsave("p_conditional_1.pdf", plot = p_conditional_1, height = 3, width = 6)

ggsave("p_conditional_2.pdf", plot = p_conditional_2, height = 3, width = 6)

# 
# fix: ensure cache path is length 1 by using a single start/end for the window

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2); ggplot2::theme_set(ggplot2::theme_bw())
  library(scales)
  library(locuszoomr)
  library(tidyr)
})

# config assumed defined: gene_symbol, flank_bp, ensdb_pkg, ld_build, ld_token, pop_codes

# read GWAS
gwas <- readr::read_tsv("print_significant_snp.txt", show_col_types = FALSE) |>
  dplyr::transmute(
    chrom = as.character(Chr),
    pos = as.numeric(bp),
    rsid = as.character(SNP),
    p = as.numeric(p),
    neglogp = -log10(p)
  ) |>
  dplyr::filter(is.finite(neglogp), !is.na(rsid), !is.na(pos))

# locus window
loc <- locuszoomr::locus(
  data   = gwas |>
    dplyr::select(chrom, pos, rsid, p),
  gene   = gene_symbol,
  flank  = flank_bp,
  ens_db = ensdb_pkg
)

# define a single region window from data
region <- gwas |>
  dplyr::filter(chrom == as.character(loc$data$chrom)) |>
  dplyr::filter(pos >= min(loc$TX$start, na.rm = TRUE),
                pos <= max(loc$TX$end,   na.rm = TRUE)) |>
  dplyr::arrange(pos)

# fallback if TX ranges are vectorised
if (nrow(region) == 0) {
  region <- gwas |>
    dplyr::filter(chrom == as.character(loc$data$chrom)) |>
    dplyr::filter(pos >= min(gwas$pos, na.rm = TRUE),
                  pos <= max(gwas$pos, na.rm = TRUE)) |>
    dplyr::arrange(pos)
}

loc_chr <- unique(region$chrom)[1]
win <- range(region$pos, na.rm = TRUE)
win_start <- as.integer(floor(win[1]))
win_end   <- as.integer(ceiling(win[2]))

# choose lead
lead_rsid <- region |>
  dplyr::arrange(p) |>
  dplyr::slice(1) |>
  dplyr::pull(rsid)

# rs list for LDmatrix
rs_list <- unique(region$rsid)
if (length(rs_list) > 300) {
  idx <- order(abs(region$pos - region$pos[match(lead_rsid, region$rsid)]))
  rs_list <- unique(region$rsid[idx])[seq_len(300)]
}

# per-population LD cache
fetch_ldmat <- function(rs_list, pop, build, token) {
  cache_dir <- "ld_cache_grid"
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  # key <- sprintf("LDmatrix_%s_%s.rds", pop, build, length(rs_list))
  key <- paste0("LDmatrix_", pop, "_", build, "_n27",#, length(rs_list),
                ".rds")
  path <- file.path(cache_dir, key)
  if (file.exists(path)) {
    readRDS(path)
  } else {
    m <- LDlinkR::LDmatrix(
      snps = rs_list,
      pop = pop,
      r2d = "r2",
      genome_build = build,
      token = token,
      file = FALSE
    )
    saveRDS(m, path)
    m
  }
}

# combined grid_df cache with single path
grid_cache_dir <- "ld_cache_grid"
if (!dir.exists(grid_cache_dir)) dir.create(grid_cache_dir, recursive = TRUE, showWarnings = FALSE)
grid_key <- sprintf("grid_%s_chr%s_%d_%d_lead_%s_%s_n%d.rds",
                    gene_symbol, loc_chr, win_start, win_end, lead_rsid, ld_build, length(rs_list))
grid_cache_path <- file.path(grid_cache_dir, grid_key)

if (file.exists(grid_cache_path)) {
  grid_df <- readRDS(grid_cache_path)
} else {
  ld_bins <- c("<0.2","0.2–0.4","0.4–0.6","0.6–0.8",">=0.8")
  allp <- lapply(pop_codes, function(pop) {
    ldmat <- try(fetch_ldmat(rs_list, pop, ld_build, ld_token), silent = TRUE)
    if (inherits(ldmat, "try-error") || is.null(ldmat)) return(NULL)
    rs_col <- dplyr::coalesce(ldmat$RS_number, ldmat[[1]])
    rownames(ldmat) <- rs_col
    colnames(ldmat)[1] <- "RS_number"
    if (!(lead_rsid %in% colnames(ldmat))) return(NULL)
    r2_vec <- suppressWarnings(as.numeric(unlist(ldmat[, lead_rsid])))
    names(r2_vec) <- rs_col
    region |>
      dplyr::left_join(
        tibble::tibble(rsid = names(r2_vec), r2 = as.numeric(r2_vec)),
        by = "rsid"
      ) |>
      dplyr::mutate(
        r2 = dplyr::if_else(rsid == lead_rsid & is.na(r2), 1, r2),
        ld_bin = cut(r2, breaks = c(-Inf,0.2,0.4,0.6,0.8,Inf), labels = ld_bins, right = FALSE),
        pop = pop
      )
  })
  grid_df <- dplyr::bind_rows(allp) |>
    dplyr::mutate(
      ld_bin = as.character(ld_bin),
      ld_bin = dplyr::if_else(is.na(ld_bin), "<0.2", ld_bin),
      ld_bin = factor(ld_bin, levels = c("<0.2","0.2–0.4","0.4–0.6","0.6–0.8",">=0.8"))
    )
  saveRDS(grid_df, grid_cache_path)
}

# combined plot coloured by population with a line across points
p_pop_r2 <- ggplot2::ggplot(
  grid_df %>% dplyr::arrange(pop, pos),
  ggplot2::aes(x = pos / 1e6, y = r2, colour = pop, group = pop)
) +
  ggplot2::geom_line(linewidth = 0.6, alpha = 0.6) +
  ggplot2::geom_point(size = 2, alpha = 0.8) +
  ggplot2::labs(
    x = paste0("position on chr", loc_chr, " (Mb)"),
    # y = expression(-log[10](italic(p))),
    title = paste0("regional association across populations relative to ", lead_rsid)
  ) +
  ggplot2::scale_colour_viridis_d(name = "population") +
  ggplot2::theme(
    legend.position = "top",
    legend.text = ggplot2::element_text(size = 8),
    legend.title = ggplot2::element_text(size = 9)
  )

p_pop_r2




# ggsave("p_pop_r2.pdf",
       # plot = p_pop_r2, height = 4, width = 8)




# standalone script: regional association across populations relative to lead, with robust caching ----
# standalone script: regional LD to lead across populations with robust caching and safe ld_bin handling
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2); ggplot2::theme_set(ggplot2::theme_bw())
  library(scales)
  library(tidyr)
  library(tibble)
  library(locuszoomr)
})

# config with safe defaults, respects pre-set variables if present
gene_symbol <- get0("gene_symbol", ifnotfound = "FAM206A")
flank_bp    <- get0("flank_bp",    ifnotfound = 1e5)
ensdb_pkg   <- get0("ensdb_pkg",   ifnotfound = "EnsDb.Hsapiens.v75")
ld_build    <- tolower(get0("ld_build", ifnotfound = "grch37"))
# ld_token    <- { t <- get0("ld_token", ifnotfound = ""); if (!nzchar(t)) Sys.getenv("LDLINK_TOKEN", ""); else t }
pop_codes   <- get0("pop_codes", ifnotfound = c(
  # "ALL",
  "AFR","AMR","EAS","EUR","SAS",
  "CEU","TSI","FIN","GBR","IBS",
  "YRI","LWK","GWD","MSL","ESN","ACB","ASW",
  "CHB","JPT","CHS","CDX","KHV",
  "MXL","PUR","CLM","PEL",
  "GIH","PJL","BEB","STU","ITU"
))

# read GWAS
gwas <- readr::read_tsv("print_significant_snp.txt", show_col_types = FALSE) %>%
  dplyr::transmute(
    chrom   = as.character(.data$Chr),
    pos     = as.numeric(.data$bp),
    rsid    = as.character(.data$SNP),
    p       = as.numeric(.data$p),
    neglogp = -log10(.data$p)
  ) %>%
  dplyr::filter(is.finite(.data$neglogp), !is.na(.data$rsid), !is.na(.data$pos))
stopifnot(nrow(gwas) > 0)

# locus window from gene and flank
loc <- locuszoomr::locus(
  data   = gwas %>% dplyr::select(.data$chrom, .data$pos, .data$rsid, .data$p),
  gene   = gene_symbol,
  flank  = flank_bp,
  ens_db = ensdb_pkg
)

# define a single region window
region <- gwas %>%
  dplyr::filter(.data$chrom == as.character(loc$data$chrom)) %>%
  dplyr::filter(
    .data$pos >= min(loc$TX$start, na.rm = TRUE),
    .data$pos <= max(loc$TX$end,   na.rm = TRUE)
  ) %>%
  dplyr::arrange(.data$pos)

if (nrow(region) == 0) {
  region <- gwas %>%
    dplyr::filter(.data$chrom == as.character(loc$data$chrom)) %>%
    dplyr::filter(
      .data$pos >= min(gwas$pos, na.rm = TRUE),
      .data$pos <= max(gwas$pos, na.rm = TRUE)
    ) %>%
    dplyr::arrange(.data$pos)
}
stopifnot(nrow(region) > 0)

# canonical keys
loc_chr   <- as.character(unique(region$chrom)[1])
win       <- range(region$pos, na.rm = TRUE)
win_start <- as.integer(floor(win[1]))
win_end   <- as.integer(ceiling(win[2]))

# lead variant
lead_rsid <- region %>%
  dplyr::arrange(.data$p) %>%
  dplyr::slice(1) %>%
  dplyr::pull(.data$rsid)

# rs list for LDmatrix, capped at 300 nearest to lead
rs_list <- unique(region$rsid)
if (length(rs_list) > 300) {
  idx_lead <- match(lead_rsid, region$rsid)
  idx <- order(abs(region$pos - region$pos[idx_lead]))
  rs_list <- unique(region$rsid[idx])[seq_len(300)]
}

# LD cache directory and helpers
ld_dir <- "ld_cache_grid"
if (!dir.exists(ld_dir)) dir.create(ld_dir, recursive = TRUE, showWarnings = FALSE)

load_ldmat_from_cache <- function(pop, build, nlen) {
  path_legacy <- file.path(ld_dir, sprintf("LDmatrix_%s_%s_n27.rds", pop, build))
  path_size   <- file.path(ld_dir, sprintf("LDmatrix_%s_%s_n%d.rds", pop, build, nlen))
  if (file.exists(path_legacy)) return(readRDS(path_legacy))
  if (file.exists(path_size))   return(readRDS(path_size))
  NULL
}

download_ldmat <- function(rsids, pop, build, token) {
  if (!requireNamespace("LDlinkR", quietly = TRUE)) return(NULL)
  if (!nzchar(token)) return(NULL)
  m <- try(LDlinkR::LDmatrix(
    snps         = rsids,
    pop          = pop,
    r2d          = "r2",
    genome_build = build,
    token        = token,
    file         = FALSE
  ), silent = TRUE)
  if (inherits(m, "try-error") || is.null(m)) return(NULL)
  saveRDS(m, file.path(ld_dir, sprintf("LDmatrix_%s_%s_n27.rds", pop, build)))
  m
}

make_reg_for_pop <- function(pop) {
  ldmat <- load_ldmat_from_cache(pop, ld_build, length(rs_list))
  if (is.null(ldmat)) ldmat <- download_ldmat(rs_list, pop, ld_build, ld_token)
  if (is.null(ldmat)) return(NULL)
  
  rs_col <- dplyr::coalesce(ldmat$RS_number, ldmat[[1]])
  rownames(ldmat) <- rs_col
  colnames(ldmat)[1] <- "RS_number"
  
  anchor <- if (lead_rsid %in% colnames(ldmat)) {
    lead_rsid
  } else {
    inter <- intersect(colnames(ldmat), region$rsid)
    if (length(inter) == 0) return(NULL)
    inter[which.min(region$p[match(inter, region$rsid)])]
  }
  
  r2_vec <- suppressWarnings(as.numeric(unlist(ldmat[, anchor])))
  names(r2_vec) <- rs_col
  
  region %>%
    dplyr::left_join(
      tibble::tibble(rsid = names(r2_vec), r2 = as.numeric(r2_vec)),
      by = "rsid"
    ) %>%
    dplyr::mutate(
      r2     = dplyr::if_else(.data$rsid == anchor & is.na(.data$r2), 1, .data$r2),
      pop    = pop,
      anchor = anchor
    )
}

# combined grid cache with single path
grid_cache_dir <- ld_dir
if (!dir.exists(grid_cache_dir)) dir.create(grid_cache_dir, recursive = TRUE, showWarnings = FALSE)
grid_key <- sprintf("grid_%s_chr%s_%d_%d_lead_%s_%s_n%d.rds",
                    gene_symbol, loc_chr, win_start, win_end, lead_rsid, ld_build, length(rs_list))
grid_cache_path <- file.path(grid_cache_dir, grid_key)

if (file.exists(grid_cache_path)) {
  grid_df <- readRDS(grid_cache_path)
} else {
  allp <- lapply(pop_codes, make_reg_for_pop)
  grid_df <- dplyr::bind_rows(allp)
  if (is.null(grid_df) || nrow(grid_df) == 0) {
    stop("No LD data available from cached or downloaded LDmatrix files for this window.")
  }
  ld_bins <- c("<0.2","0.2–0.4","0.4–0.6","0.6–0.8",">=0.8")
  grid_df <- grid_df %>%
    dplyr::mutate(
      ld_bin_chr = cut(.data$r2, breaks = c(-Inf,0.2,0.4,0.6,0.8,Inf), labels = ld_bins, right = FALSE),
      ld_bin_chr = as.character(.data$ld_bin_chr),
      ld_bin_chr = dplyr::if_else(is.na(.data$ld_bin_chr), "<0.2", .data$ld_bin_chr),
      ld_bin     = factor(.data$ld_bin_chr, levels = ld_bins, ordered = TRUE)
    ) %>%
    dplyr::select(-ld_bin_chr)
  saveRDS(grid_df, grid_cache_path)
}


# draw lines that follow the same points per population
grid_df_line <- grid_df %>%
  dplyr::filter(is.finite(r2)) %>%
  dplyr::group_by(pop) %>%
  dplyr::arrange(pos, .by_group = TRUE) %>%
  dplyr::ungroup()

p_pop_r2 <- ggplot() +
  ggplot2::geom_path(
    data = grid_df_line,
    ggplot2::aes(x = pos / 1e6, y = r2, colour = pop, group = pop),
    linewidth = 0.6, alpha = 0.6, na.rm = TRUE
  ) +
  ggplot2::geom_point(
    data = grid_df,
    ggplot2::aes(x = pos / 1e6, y = r2, colour = pop),
    size = 2, alpha = 0.8, na.rm = TRUE
  ) +
  ggplot2::labs(
    x = paste0("position on chr", loc_chr, " (Mb)"),
    y = expression(LD~r^2~to~anchor),
    title = paste0("regional LD across populations")
    # relative to lead ", lead_rsid)
  ) +
  ggplot2::scale_colour_viridis_d(name = "population", option = "C") +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::theme(
    legend.position = "top",
    legend.text = ggplot2::element_text(size = 7),
    legend.title = ggplot2::element_text(size = 9)
  )

print(p_pop_r2)
# ggsave("p_pop_r2.pdf", plot = p_pop_r2, height = 4, width = 8)




# add square outline spanning full y and x = 111,802,731 to 111,837,776
xmin_bp <- min(111737776, 111832731)
xmax_bp <- max(111737776, 111832731)
ymin_r2 <- min(grid_df$r2, na.rm = TRUE)
ymax_r2 <- max(grid_df$r2, na.rm = TRUE)

p_pop_r2 <- p_pop_r2 +
  ggplot2::annotate(
    "rect",
    xmin = xmin_bp / 1e6, xmax = xmax_bp / 1e6,
    ymin = ymin_r2,       ymax = ymax_r2,
    fill = NA, colour = "black", linewidth = 0.7
  ) 

print(p_pop_r2)




xmin_bp <- min(111660776, 111730776)
xmax_bp <- max(111660776, 111730776)
ymin_r2 <- min(grid_df$r2, na.rm = TRUE)
ymax_r2 <- max(grid_df$r2, na.rm = TRUE)

p_pop_r2 <- p_pop_r2 +
  ggplot2::annotate(
    "rect",
    xmin = xmin_bp / 1e6, xmax = xmax_bp / 1e6,
    ymin = ymin_r2,       ymax = ymax_r2,
    fill = NA, colour = "black", linewidth = 0.7
  ) 

p_pop_r2 <- p_pop_r2 +
ggplot2::theme(
  legend.position = "right",
  legend.text = ggplot2::element_text(size = 7),
  legend.title = ggplot2::element_text(size = 9)
) +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 4)) 



p_pop_r2 <- p_pop_r2 + xlim(111.55, 111.95)


print(p_pop_r2)
ggsave("p_pop_r2_box.pdf", plot = p_pop_r2, height = 2.5, width = 8)





# Susie coloc ----


library(dplyr)
library(readr)
library(LDlinkR)
library(susieR)
library(Matrix)

# read GWAS
gwas <- read_tsv("print_significant_snp.txt", show_col_types = FALSE) %>%
  transmute(
    rsid = as.character(SNP),
    pos  = as.numeric(bp),
    beta = as.numeric(b),
    se   = as.numeric(se),
    p    = as.numeric(p)
  )

# choose lead SNP
lead <- gwas %>% arrange(p) %>% slice(1) %>% pull(rsid)

# list rsIDs
rs_list <- unique(gwas$rsid)

# request reference LD
ldmat <- LDmatrix(
  snps = rs_list,
  pop = "EUR",
  r2d = "r2",
  token = ld_token,
  genome_build = ld_build,
  file = FALSE
)

# tidy LD into numeric matrix
rs_col <- coalesce(ldmat$RS_number, ldmat[[1]])
rownames(ldmat) <- rs_col
colnames(ldmat)[1] <- "RS_number"
R <- as.matrix(ldmat[, -1, drop = FALSE])
rownames(R) <- rs_col
colnames(R) <- colnames(R)

# match GWAS to LD
g <- gwas %>% filter(rsid %in% rownames(R))
R <- R[g$rsid, g$rsid, drop = FALSE]
R <- nearPD(R, corr = TRUE)$mat %>% as.matrix()

# z and N
z <- g$beta / g$se
N <- 510 + 994

# run SuSiE
fit <- susie_rss(
  z = z,
  R = R,
  n = N,
  L = 3,
  max_iter = 200,
  coverage = 0.95,
  verbose = FALSE
)

# extract credible sets
cs <- susie_get_cs(fit)

# output results
list(
  credible_sets = cs$cs,
  pip = fit$pip,
  lead_variants = sapply(cs$cs, function(idx) g$rsid[idx[which.max(fit$pip[idx])]])
)


## GWAS z vs GTEx eQTL beta ----


library(readxl)
library(dplyr)
library(ggplot2)
library(readxl)
library(dplyr)
library(ggplot2)

# read sheets, skipping the first row (table title)
eqtl <- read_excel(
  "lawless2025spss_gwas_tables.xlsx",
  sheet = "S5.eQTL",
  skip = 1
)

sqtl <- read_excel(
  "lawless2025spss_gwas_tables.xlsx",
  sheet = "S6.sQTL",
  skip = 1
)

names(sqtl)

# clean up column names
eqtl_clean <- eqtl %>%
  transmute(
    rsid = `SNP Id`,
    gene = `Gene Symbol`,
    tissue = Tissue,
    nes = NES,
    p = `P-Value`,
    type = "eQTL"
  )

sqtl_clean <- sqtl %>%
  transmute(
    rsid = `SNP Id`,
    gene = `Gene Symbol`,
    tissue = Tissue,
    nes = NES,
    p = `P-Value`,
    type = "sQTL"
  )

# combine both
gtex_all <- bind_rows(eqtl_clean, sqtl_clean)

# keep only GWAS locus SNPs
gtex_filtered <- gtex_all %>% filter(rsid %in% g$rsid)

### Co-localisation plot option 1: GWAS z vs GTEx NES ----

coloc_df <- g %>%
  mutate(z = beta / se) %>%
  inner_join(gtex_filtered, by = "rsid")

ggplot(coloc_df, aes(x = z, y = nes, colour = type, shape = type)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  theme_bw() +
  labs(
    x = "GWAS z score (sepsis susceptibility)",
    y = "GTEx NES (effect on expression or splicing)",
    colour = "",
    shape = ""
  )


### Co-localisation plot option 2: PIP vs genomic position ----

pip_df <- g %>%
  mutate(
    pip = fit$pip[match(rsid, names(fit$pip))]
  )

plot_pip <- ggplot(pip_df, aes(x = pos, y = pip)) +
  geom_point(size = 2) +
  geom_point(
    data = pip_df %>% filter(pip == max(pip)),
    colour = "red", size = 3
  ) +
  theme_bw() +
  labs(
    x = "Genomic position (bp)",
    y = "Posterior\ninclusion probability"
  )

plot_pip

ggsave("p_plot_pip.pdf", plot_pip, width = 8, height = 2.5)

print("This produces a clean plot showing a single red point at rs28361152 with PIP = 1.0 and all other variants at zero.")
### across tiussue ----

library(readxl)
library(dplyr)
library(ggplot2)

eqtl <- read_excel(
  "lawless2025spss_gwas_tables.xlsx",
  sheet = "S5.eQTL",
  skip = 1
)

sqtl <- read_excel(
  "lawless2025spss_gwas_tables.xlsx",
  sheet = "S6.sQTL",
  skip = 1
)

eqtl_clean <- eqtl %>%
  transmute(
    rsid = `SNP Id`,
    gene = `Gene Symbol`,
    tissue = Tissue,
    nes = NES,
    type = "eQTL"
  )

sqtl_clean <- sqtl %>%
  transmute(
    rsid = `SNP Id`,
    gene = `Gene Symbol`,
    tissue = Tissue,
    nes = NES,
    type = "sQTL"
  )

gtex_all <- bind_rows(eqtl_clean, sqtl_clean)

gtex_rs <- gtex_all %>% filter(rsid == "rs28361152")

gtex_rs <- unique(gtex_rs)


genes_left  <- unique(df$gene)[1:2]
genes_right <- unique(df$gene)[3:5]

df_left  <- df  %>% filter(gene %in% genes_left)
df_right <- df %>% filter(gene %in% genes_right)

labels_left  <- labels_df %>% filter(gene %in% genes_left)
labels_right <- labels_df %>% filter(gene %in% genes_right)

p_left <- ggplot(
  df_left,
  aes(
    x = reorder(tissue, nes),
    y = nes,
    fill = type,
    colour = type
  )
) +
  geom_col(width = 0.6, linewidth = 0.3) +
  geom_label(
    data = labels_left,
    aes(x = label_x, y = label_y, label = gene),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1,
    fill = scales::alpha("white", 0.5),
    colour = "black",
    label.size = 0.4,
    label.padding = unit(0.2, "lines"),
    size = 3
  ) +
  scale_fill_manual(values = c("eQTL" = "#1f78b4", "sQTL" = "#e31a1c")) +
  scale_colour_manual(values = c("eQTL" = "black", "sQTL" = "black")) +
  coord_flip(clip = "off") +
  facet_grid(
    rows = vars(gene),
    scales = "free_y",
    space = "free_y"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 7),
    strip.text = element_blank(),
    panel.spacing = unit(0.5, "lines")
  ) +  labs(
    x = "Tissue",
    y = "GTEx NES for rs28361152"
  )


p_left

p_right <- ggplot(
  df_right,
  aes(
    x = reorder(tissue, nes),
    y = nes,
    fill = type,
    colour = type
  )
) +
  geom_col(width = 0.6, linewidth = 0.3) +
  geom_label(
    data = labels_right,
    aes(x = label_x, y = label_y, label = gene),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1,
    fill = scales::alpha("white", 0.5),
    colour = "black",
    label.size = 0.4,
    label.padding = unit(0.2, "lines"),
    size = 3
  ) +
  scale_fill_manual(values = c("eQTL" = "#1f78b4", "sQTL" = "#e31a1c")) +
  scale_colour_manual(values = c("eQTL" = "black", "sQTL" = "black")) +
  coord_flip(clip = "off") +
  facet_grid(
    rows = vars(gene),
    scales = "free_y",
    space = "free_y"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 7),
    strip.text = element_blank(),
    panel.spacing = unit(0.5, "lines")
  ) +  labs(
    x = "Tissue",
    y = "GTEx NES for rs28361152"
  )

p_right

library(patchwork)

p_final <- p_left + p_right + plot_layout(ncol = 2) 
p_final
ggsave("p_gtex_NES_split2col.pdf", p_final, width = 8, height = 5)

