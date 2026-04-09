rm(list=ls())

# an analysis investigating the association between the abundance of any specific genus at baseline and later irae

# Define your list of packages
packages <- c(
  "RColorBrewer", "ggplot2", "grid", "MASS", "pheatmap", "viridis", "vegan", "ape",
  "Rtsne", "umap", "GGally", "gtools", "tidyr", "readxl", "ggpubr",
  "rstatix", "tidyverse", "survival", "survminer", "tibble", "stringr", "dplyr", "gridExtra"
)

# Install any packages that aren't already installed
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

# Load all packages
lapply(packages, library, character.only = TRUE)

# set the working directory as the location of the script ---
setwd("S:/MolecularOncology/Users/Michael C/neoadjuvant microbiome/NGS results/microbiome_R_directory/indir/")
genera <- read.table("jhmi-neoadj-mbiom-2024.genus_baseline.txt", sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)

genera$any_irae <- as.factor(genera$any_irae)
genera$high_grade_irae <- as.factor(genera$high_grade_irae)
genera$rad_response_preop <- as.factor(genera$rad_response_preop)

# Create a results data frame
iraeresults <- data.frame(
  Genus = colnames(genera)[30:ncol(genera)],
  p_value = NA,
  stringsAsFactors = FALSE
)

# Assume results table has a column 'Genus' with full names
iraeresults$Genus <- sub(".*g_", "", iraeresults$Genus)

# analysis for irae

# Loop through each genus
for(i in 30:ncol(genera)){
  genus_values <- genera[, i]
  iraeresults$p_value[i - 29] <- wilcox.test(genus_values ~ genera$any_irae)$p.value
}


# View significant genera
iraesignificant_genera <- iraeresults[iraeresults$p_value < 0.05, ]
iraesignificant_genera

# Remove NAs, keep only significant genera, and sort
iraesignificant_genera_sorted <- iraeresults[!is.na(iraeresults$p_value) & iraeresults$p_value < 0.05, ]  # filter
iraesignificant_genera_sorted <- iraesignificant_genera_sorted[order(iraesignificant_genera_sorted$p_value), ]  # sort

# View the final table
iraesignificant_genera_sorted


# analysis for hgirae

# Create a results data frame
hgiraeresults <- data.frame(
  Genus = colnames(genera)[30:ncol(genera)],
  p_value = NA,
  stringsAsFactors = FALSE
)

# Assume results table has a column 'Genus' with full names
hgiraeresults$Genus <- sub(".*g_", "", hgiraeresults$Genus)

# Loop through each genus
for(i in 30:ncol(genera)){
  genus_values <- genera[, i]
  hgiraeresults$p_value[i - 29] <- wilcox.test(genus_values ~ genera$high_grade_irae)$p.value
}


# View significant genera
hgiraeresultssignificant_genera <- hgiraeresults[hgiraeresults$p_value < 0.05, ]
hgiraeresultssignificant_genera

# Remove NAs, keep only significant genera, and sort
hgiraeresultssignificant_genera_sorted <- hgiraeresults[!is.na(hgiraeresults$p_value) & hgiraeresults$p_value < 0.05, ]  # filter
hgiraeresultssignificant_genera_sorted <- hgiraeresultssignificant_genera_sorted[order(hgiraeresultssignificant_genera_sorted$p_value), ]  # sort

# View the final table
hgiraeresultssignificant_genera_sorted


# Find max number of rows
max_rows <- max(nrow(iraesignificant_genera_sorted), nrow(hgiraeresultssignificant_genera_sorted))

# Pad MPR table if needed
if(nrow(iraesignificant_genera_sorted) < max_rows){
  iraesignificant_genera_sorted[(nrow(iraesignificant_genera_sorted)+1):max_rows, ] <- NA
}

# Pad RR table if needed
if(nrow(hgiraeresultssignificant_genera_sorted) < max_rows){
  hgiraeresultssignificant_genera_sorted[(nrow(hgiraeresultssignificant_genera_sorted)+1):max_rows, ] <- NA
}

combined_tables <- cbind(
  iraesignificant_genera_sorted,
  hgiraeresultssignificant_genera_sorted
)

# Optional: rename columns to indicate outcome
colnames(combined_tables) <- c(
  paste0("irae_", colnames(iraesignificant_genera_sorted)),
  paste0("hgirae_", colnames(hgiraeresultssignificant_genera_sorted))
)

# View the side-by-side table
combined_tables

# Load necessary packages
library(gridExtra)

# Save to PDF

# Create the table grob
table_grob <- tableGrob(combined_tables)

# Create the title grob
title <- textGrob(
  "Baseline gut microbiome genera significantly associated with immune related adverse events",
  gp = gpar(fontsize = 14, fontface = "bold")
)

# Save to PDF
pdf("combined_tables_with_title.pdf", width = 11, height = 8.5)
grid.arrange(title, table_grob, ncol = 1, heights = c(0.1, 0.9))
dev.off()
