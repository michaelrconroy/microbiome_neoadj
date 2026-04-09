rm(list=ls())
setwd("S:/MolecularOncology/Users/Michael C/neoadjuvant microbiome/NGS results/microbiome_R_directory")
analysisdir = "./analysis/heatmaps/"
getwd()

# Define your list of packages
packages <- c(
  "RColorBrewer", "ggplot2", "grid", "MASS", "pheatmap", "viridis", "vegan", "ape",
  "Rtsne", "umap", "GGally", "gtools", "tidyr", "dplyr", "readxl", "ggpubr",
  "rstatix", "tidyverse", "survival", "survminer", "broom", "tibble", "stringr"
)

# Install any packages that aren't already installed
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

# Load all packages
lapply(packages, library, character.only = TRUE)

df_abund <- read_excel("jhmi-neoadj-mbiom-2024.genus_baseline_10d.xlsx")

# Load required libraries
library(dplyr)
library(tibble)

# 1) Extract original genus columns (e.g., columns 29 to 259)
genus_cols <- names(df_abund)[29:259]

# 2) Extract just genus names (everything after "g_")
genus_short <- str_extract(genus_cols, "(?<=g_)[^\\.]+")  

# 3) Replace any NA (from missing "g_") with original name to be safe
genus_short[is.na(genus_short)] <- genus_cols[is.na(genus_short)]

# 4) Remove unwanted genus names from the short list
keep_idx <- !(genus_short %in% c("unassigned", "Incertae_sedis"))
genus_short_filtered <- genus_short[keep_idx]
genus_cols_filtered <- genus_cols[keep_idx]

# 5) Rename the columns in df_abund (only the filtered ones)
names(df_abund)[29:259][keep_idx] <- genus_short_filtered

# 6) Drop the unwanted columns from df_abund
df_abund <- df_abund %>%
  select(-all_of(genus_cols[!keep_idx]))

# 7) Update genus_cols to the filtered & truncated names for your downstream code
genus_cols <- genus_short_filtered


# Perform Wilcoxon rank-sum test for each genus
results <- lapply(genus_cols, function(genus) {
  group1 <- df_abund %>% filter(MPR == 0) %>% pull(!!sym(genus))
  group2 <- df_abund %>% filter(MPR == 1) %>% pull(!!sym(genus))
  
  test <- wilcox.test(group1, group2, exact = FALSE)
  
  tibble(
    genus = genus,
    median_group1 = median(group1, na.rm = TRUE),
    median_group2 = median(group2, na.rm = TRUE),
    p_value = test$p.value,
    statistic = test$statistic
  )
})

# Combine into one dataframe and sort by p-value
mw_results <- bind_rows(results) %>%
  arrange(p_value)

mw_results


sig_genera <- mw_results %>%
  filter(p_value < 0.05) %>%
  pull(genus)

# 2. Long format for plotting
df_long <- df_abund %>%
  select(MPR, all_of(sig_genera)) %>%
  tidyr::pivot_longer(
    cols = -MPR,
    names_to = "Genus",
    values_to = "Abundance"
  ) %>%
  mutate(MPR = factor(MPR, levels = c(0, 1), labels = c("non-MPR", "MPR")))

# 3. Calculate p-value label positions for each genus
label_positions <- df_long %>%
  group_by(Genus) %>%
  summarise(label_y = max(Abundance, na.rm = TRUE) * 1.1, .groups = "drop")

# 4. Merge label positions into the plotting data
df_long <- df_long %>%
  left_join(label_positions, by = "Genus")

# 5. Plot
p <- ggplot(df_long, aes(x = MPR, y = Abundance, fill = MPR)) +
  geom_boxplot(outlier.shape = NA) +
  # Points with Abundance == 0, no jitter
  geom_point(data = subset(df_long, Abundance == 0), alpha = 0.6, size = 1)+
  # Points with Abundance != 0, jitter applied
  geom_jitter(data = subset(df_long, Abundance != 0), 
              width = 0.2, alpha = 0.6, size = 1) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    aes(label.y = label_y), size = 3
  ) +
  facet_wrap(~ Genus, scales = "free_y") +
  labs(x = "MPR status", y = "Relative abundance", fill = "MPR") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"))

ggsave("boxplots_significant_genera_MPR.pdf", plot = p, width = 12, height = 8)
