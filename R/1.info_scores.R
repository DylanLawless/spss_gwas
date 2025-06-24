library(ggplot2)
library(dplyr)
library(readr)

# Load imputation info score data
info_sample <- read_table2(
  file = "./SHCS.QC_chr9.pos96347718-101347219.impute2_snptest_score",
  na = c("", "NA")
)

# Preview
head(info_sample, 100)

# Histogram of INFO scores
ggplot(info_sample, aes(x = as.numeric(info))) +
  geom_histogram(binwidth = 0.01, fill = "#3366CC", colour = "white") +
  labs(
    title = "Distribution of Imputation INFO Scores",
    x = "INFO Score",
    y = "Count"
  ) +
  theme_minimal()

# Cumulative frequency plot
info_sample %>%
  filter(info >= 0) %>%
  count(info) %>%
  arrange(info) %>%
  mutate(cum_n = cumsum(n)) %>%
  ggplot(aes(x = info, y = cum_n)) +
  geom_line(colour = "#CC0000") +
  geom_point(size = 1.5) +
  labs(
    title = "Cumulative Frequency of INFO Scores",
    x = "INFO Score",
    y = "Cumulative Variant Count"
  ) +
  theme_minimal()

