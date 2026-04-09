rm(list=ls())

# Define your list of packages
packages <- c(
  "RColorBrewer", "ggplot2", "grid", "MASS", "pheatmap", "viridis", "vegan", "ape",
  "Rtsne", "umap", "GGally", "gtools", "tidyr", "readxl", "ggpubr",
  "rstatix", "tidyverse", "survival", "survminer", "tibble", "stringr", "dplyr"
)

# Install any packages that aren't already installed
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

# Load all packages
lapply(packages, library, character.only = TRUE)

## identifying if differences in alpha diversity between pre and post specimens

# set the working directory as the location of the script ---
setwd("S:/MolecularOncology/Users/Michael C/neoadjuvant microbiome/NGS results/microbiome_R_directory/indir/")
A <- read.table("jhmi-neoadj-mbiom-2024.alpha_pre.txt", sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)
B <- read.table("jhmi-neoadj-mbiom-2024.alpha_post.txt", sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)
A$timepoint <- "before"
B$timepoint <- "after"
combined <- rbind(A, B)
combined <- combined %>%
  dplyr::select(SampleID, timepoint, everything())


# Identify the alpha diversity columns (adjust column indices or names)
alpha_cols <- c("simpson_reciprocal", "shannon", "chao1", "observed_species")

# Initialize results dataframe with extra columns for medians
results <- data.frame(
  metric = alpha_cols,
  before_median = NA,
  after_median = NA,
  test_statistic = NA,
  p_value = NA
)

# Loop through each metric
for (i in seq_along(alpha_cols)) {
  
  metric <- alpha_cols[i]
  
  # Convert to wide format so each row is one patient
  df_wide <- combined %>%
    dplyr::select(Record_ID, timepoint, all_of(metric)) %>%
    tidyr::pivot_wider(names_from = timepoint, values_from = all_of(metric))
  
  # Perform paired Wilcoxon signed-rank test
  test <- wilcox.test(df_wide$before, df_wide$after, paired = TRUE)
  
  # Compute medians for each timepoint
  before_median <- median(df_wide$before, na.rm = TRUE)
  after_median  <- median(df_wide$after, na.rm = TRUE)
  
  # Store results
  results$before_median[i] <- before_median
  results$after_median[i]  <- after_median
  results$test_statistic[i] <- test$statistic
  results$p_value[i] <- test$p.value
}

results

library(openxlsx)

# Round and clean up your results first
simple_results <- results %>%
  dplyr::mutate(
    before_median = round(before_median, 2),
    after_median  = round(after_median, 2),
    test_statistic = round(test_statistic, 2),
    p_value = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  )

# Save as Excel file
write.xlsx(simple_results, file = "alpha_diversity_results.xlsx", rowNames = FALSE)