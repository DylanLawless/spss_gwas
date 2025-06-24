require(ggplot2)
require(dplyr)
require(gridExtra)

#=========================================================
# Eigen value plot
#=========================================================
val <- read.table("./combined_grm/cohort.filt.pca.eigenval", header=F)
val <- val %>% mutate(V2 = rownames(val))

val %>%
  ggplot(aes(x = as.numeric(V2), y = V1))+
  geom_point()

#=========================================================
# Eigen vectors
#=========================================================
spss.fam <- read.table("./combined_grm/post_qc_sample_list/SPSS.fam", header=F)
shcs.fam <- read.table("./combined_grm/post_qc_sample_list/SHCS.fam", header=F)
vec <- read.table("./combined_grm/cohort.filt.pca.eigenvec", header=F)

fams <- rbind(shcs.fam, spss.fam)
colnames(fams) <- c("V1", "V2", "x1", "x2", "sex", "Group")

merged <- merge(x=vec, y=fams, by=c("V1","V2"), all=F)
rm(vec, fams, spss.fam, shcs.fam)

sepsis_cohort <- read.csv("./combined_grm/spss_gwas_episode.csv", header=T, sep = ",")
colnames(sepsis_cohort)[colnames(sepsis_cohort) == "sample.id"] <- "V1"
merged_2 <- merge(x=merged, y=sepsis_cohort, by="V1", all=F)

# PCA_cases_ethnicity.pdf 5x6
merged_2 %>%
  filter(Group == 2) %>%
  ggplot(aes(x=V3, y=V4)) +
  geom_point(aes(fill = ethnicity), shape=21, alpha=1) +
  labs(x = "PC1", y = "PC2") +
  theme(text = element_text(face="bold"), legend.position="right", panel.background = element_rect("#F7F7F7"))

# 2b.PCA_case_control_EU.pdf 5x6
merged %>%
  ggplot(aes(x=V3, y=V4)) +
  geom_point(aes(fill = as.character(Group)), shape=21, alpha=1) +
  labs(x = "PC1", y = "PC2") +
  theme(text = element_text(face="bold"), legend.position="right", panel.background = element_rect("#F7F7F7")) +
  scale_fill_discrete(name = "Group", labels = c("Control", "Case")) +
  geom_vline(xintercept = -0.0055) +
  geom_hline(yintercept = 0.0)

merged_EU <- merged %>%
  filter(V3 < -0.0055, V4 > 0.0, V7 < 0.059, V13 < 0.025, V19 < 0.035, V20 < 0.013)

merged_non_EU <- merged %>%
  filter(V3 > -0.0055 |
         V4 < 0.0 |
         V7 > 0.059 |
         V13 > 0.025 |
         V19 > 0.035 |
         V20 > 0.013)

merged_non_EU %>% group_by(Group) %>% tally()
merged_EU %>% group_by(Group) %>% tally()
merged %>% group_by(Group) %>% tally()

# 2.PCA_case_control_EU.pdf 5x6
merged_EU %>%
  ggplot(aes(x=V3, y=V4)) +
  geom_point(aes(fill = as.character(Group)), shape=21, alpha=1) +
  labs(x = "PC1", y = "PC2") +
  theme(text = element_text(face="bold"), legend.position="right", panel.background = element_rect("#F7F7F7")) +
  scale_fill_discrete(name = "Group", labels = c("Control", "Case"))

merged_non_EU %>%
  select(V1, V2) %>%
  write.table(file='merged_non_EU.csv', sep=",", quote=FALSE, row.names=FALSE)

# 3.PCA_case_control_EU_20pc.pdf 12x16
pca_plots <- lapply(1:19, function(i) {
  merged_EU %>%
    ggplot(aes_string(x = paste0("V", i + 2), y = paste0("V", i + 3), colour = "as.character(Group)")) +
    geom_point(alpha = 0.5) +
    labs(x = paste0("PC", i), y = paste0("PC", i + 1)) +
    theme(text = element_text(face="bold"), legend.position = if (i == 19) "right" else "none",
          panel.background = element_rect("#F7F7F7")) +
    if (i == 19) scale_colour_discrete(name = "Group", labels = c("Control", "Case")) else NULL
})

grid.arrange(grobs = pca_plots, ncol = 4, top = "Cohort Principal component analysis", left = "")

