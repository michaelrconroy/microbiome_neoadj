rm(list=ls())

# an analysis comparing all 'before' specimens to all 'after' specimens (not necessarily paired) to see if there are any genera with significantly different abundances.

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

## identifying which genera have statistically significant in abundance between the 'pre' and 'post' specimens

# set the working directory as the location of the script ---
setwd("S:/MolecularOncology/Users/Michael C/neoadjuvant microbiome/NGS results/microbiome_R_directory/indir/")
genera <- read.table("jhmi-neoadj-mbiom-2024.genus.txt", sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)

# Create a new column for IO group
genera$IO_group <- ifelse(genera$Days_from_IO_start <= 2, "Before_IO", "After_IO")

# Assuming columns 30 to end are genus abundances
genus_data <- genera[, 30:ncol(genera)]

# Initialize a results data frame
results <- data.frame(
  Genus = colnames(genus_data),
  p_value = NA,
  stringsAsFactors = FALSE
)

# Loop through each genus
for(i in 1:ncol(genus_data)){
  genus_values <- genus_data[, i]
  results$p_value[i] <- wilcox.test(
    genus_values ~ genera$IO_group
  )$p.value
}


# View significant results
significant_genera <- results[results$padj < 0.05, ]
significant_genera

