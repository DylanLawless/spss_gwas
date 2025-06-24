require(ggplot2)
require(dplyr)
require(gridExtra)

#=========================================================
# After removing outlier groups and rerun PCA
#=========================================================

#=========================================================
# Eigen vectors
#=========================================================
spss.fam <- read.table("./combined_grm/post_qc_sample_list/SPSS.fam", header=F)
shcs.fam <- read.table("./combined_grm/post_qc_sample_list/SHCS.fam", header=F)
vec <- read.table("./combined_grm/cohort.filt.group.pca.eigenvec", header=F)

fams <- rbind(shcs.fam, spss.fam)
colnames(fams) <- c("V1", "V2", "x1", "x2", "sex", "Group")

merged_grouped <- merge(x=vec, y=fams, by=c("V1","V2"), all=F)
rm(vec, fams, spss.fam, shcs.fam)

# PCA plots for PC1–PC20
pca_plots <- lapply(1:19, function(i) {
  p <- merged_grouped %>%
    ggplot(aes_string(x = paste0("V", i + 2), y = paste0("V", i + 3), colour = "as.character(Group)")) +
    geom_point(alpha = 0.5) +
    labs(x = paste0("PC", i), y = paste0("PC", i + 1)) +
    theme(text = element_text(face="bold"),
          legend.position = if (i == 19) "right" else "none",
          panel.background = element_rect("#F7F7F7"))
  if (i == 19) {
    p <- p + scale_colour_discrete(name = "Group", labels = c("Control", "Case"))
  }
  p
})

grid.arrange(grobs = pca_plots, ncol=4, top = "Cohort Principal component analysis", left = "")
# 4.PCA_case_control_EU_20pc_after_prune.pdf 12x16

