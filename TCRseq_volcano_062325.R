rm(list=ls())

library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations
library(readxl)

# Set input path
path <- "S:/MolecularOncology/Users/Michael C/neoadjuvant microbiome/TCRseq/"
setwd(path)

# Import TCR results
df <- read_xlsx(paste0(path, 'TCR_genus_association.xlsx'),sheet=2)
df$spearman_est <- as.numeric(df$spearman_est)

#format data

df$sp_pval <- as.numeric(df$sp_pval)
df <- df %>% filter(!is.na(sp_pval))

# Add a column to the data frame to specify if they are UP- or DOWN- regulated 

df$diffexpressed <- NA
# if Spearman > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$spearman_est > 0.6 & df$sp_pval < 0.05] <- "Positive correlation, p <0.05"
# if Spearman < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$spearman_est < -0.6 & df$sp_pval < 0.05] <- "Negative correlation, p <0.05"

df$diffexpressed <- as.character(df$diffexpressed)  # Convert to character to modify
df$diffexpressed[is.na(df$diffexpressed)] <- "no significant correlation"
df$diffexpressed <- factor(df$diffexpressed, levels = c(
  "Negative correlation, p <0.05",
  "Positive correlation, p <0.05",
  "no significant correlation"
))



df$feature <- as.character(df$feature)  # Ensure character
df$spearman_est <- as.numeric(df$spearman_est)  # Ensure numeric

df$feature <- trimws(tolower(df$feature))  # Standardize formatting

# Extract the top 10 highest and bottom 10 lowest spearman_est features
top10_features <- df %>%
  arrange(desc(spearman_est)) %>%
  head(10) %>%
  pull(feature)

bottom10_features <- df %>%
  arrange(spearman_est) %>%
  head(10) %>%
  pull(feature)

# Combine both top and bottom features
selected_features <- c(top10_features, bottom10_features)

selected_features <- trimws(tolower(selected_features))  # Standardize names

# Label only selected features
df$delabel <- ifelse(df$feature %in% selected_features, df$feature, NA_character_)

# Extract only the genus name
df$genus <- sub(".*g_", "", df$delabel)
df$family <- str_extract(df$delabel, "(?<=\\.f_)[^\\.]+")


#volcano plot

theme_set(theme_classic(base_size = 12) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5),
            ))

p <- ggplot(data = df, aes(x = spearman_est, y = -log10(sp_pval), col=diffexpressed, label=genus)) + scale_y_continuous(breaks = c(1.0,2.0,3.0,4.0)) + geom_vline(xintercept = c(-0.6, 0.6), col = "red", linetype = 'dashed') +
  geom_hline(yintercept = 0.05, col = "red", linetype = 'dashed') + 
  geom_point() + labs(x = "Spearman Estimate of correlation between MPF and genus", y = "-log10 p value of Spearman Estimate") +   scale_color_manual(name = "Genera with  correlation",values = c("#8B0000", "#00008B","#808080"
  ), 
                                                                                                                                               labels = c("Decreased abundance, p <0.05", "Increased abundance, p <0.05", "No significant correlation")) +
  geom_text_repel(data = df %>% filter(feature %in% top10_features),
                  aes(label = genus),
                  max.overlaps = Inf, 
                  nudge_x = 2,       # Nudges top 10 to the right
                  direction = "y",   # Keeps vertical alignment
                  hjust = 0,                  force = 5,           # Increase repulsion force
                  box.padding = 0.5,   # More spacing between labels
                  point.padding = 0.3) +       # Align text to the left of the callout line
  
  geom_text_repel(data = df %>% filter(feature %in% bottom10_features),
                  aes(label = genus),
                  max.overlaps = Inf, 
                  nudge_x = -2,      # Nudges bottom 10 to the left
                  direction = "y",   # Keeps vertical alignment
                  hjust = 1, 
                  force = 5,           # Increase repulsion force
                  box.padding = 0.5,   # More spacing between labels
                  point.padding = 0.3)  # More spacing from points# Align text to the right of the callout line

ggsave("volcano_genus_0725.pdf", plot = p, width = 8, height = 6, dpi = 300)


#volcano plot by family

volcano <- ggplot(data = df, aes(x = spearman_est, y = -log10(sp_pval), col=diffexpressed, label=family)) + 
  geom_vline(xintercept= c(-0.6, 0.6), col = "red", linetype = 'dashed') +
  geom_hline(yintercept = 1.30, col = "red", linetype = 'dashed') + 
  geom_point() + labs(x = "Spearman Estimate of correlation between MPF and family", y = "-log10 of p value of Spearman Estimate",  title = "Significant correlations between abundances of \n family and clonality of T cell infiltrate") +   
  scale_color_manual(name = "Families with  correlation", values = c("#8B0000", "#F08080", "#ADD8E6", "#00008B"),                  labels = c("Decreased abundance, p <0.05", "Decreased abundance, Not significant", "Increased abundance, Not significant", "Increased abundance, p <0.05")) + 
  coord_cartesian(ylim = c(0, 5), xlim = c(-3, 3)) +
  geom_text_repel(data = df %>% filter(feature %in% top10_features),
                  aes(label = family),
                  max.overlaps = Inf, 
                  key_glyph = draw_key_point,
                  nudge_x = 2,       # Nudges top 10 to the right
                  direction = "y",   # Keeps vertical alignment
                  hjust = 0,                  
                  force = 12,           # Increase repulsion force
                  box.padding = 0.8,   # More spacing between labels
                  point.padding = 0.3) +       # Align text to the left of the callout line
  
  geom_text_repel(data = df %>% filter(feature %in% bottom10_features),
                  aes(label = family),
                  max.overlaps = Inf,
                  key_glyph = draw_key_point,
                  nudge_x = -2,      # Nudges bottom 10 to the left
                  direction = "y",   # Keeps vertical alignment
                  hjust = 0, 
                  force = 12,           # Increase repulsion force
                  box.padding = 0.8,   # More spacing between labels
                  point.padding = 0.3)+  # More spacing from points# Align text to the right of the callout line
  theme(
    plot.margin = margin(20, 20, 20, 50),  # Adjust the margins (top, right, bottom, left)
    plot.title = element_text(size = rel(1.1),face="bold",hjust = 0.5)  # Make title font smaller and center it
  )

ggsave("volcano_062325.png", plot = volcano, width = 8, height = 6, dpi = 300)

#volcano plot by genus

volcano <- ggplot(data = df, aes(x = spearman_est, y = -log10(sp_pval), col=diffexpressed, label=genus)) + 
  geom_vline(xintercept= c(-0.6, 0.6), col = "red", linetype = 'dashed') +
  geom_hline(yintercept = 1.30, col = "red", linetype = 'dashed') + 
  geom_point() + labs(x = "Spearman Estimate of correlation between MPF and genus", y = "-log10 of p value of Spearman Estimate",  title = "Significant correlations between abundances of \n genus and clonality of T cell infiltrate") +   
  scale_color_manual(name = "Genera with significant correlation", values = c("#00AFBB", "grey", "#006400"),                  labels = c("Decreased abundance", "Not significant", "Increased abundance")) + 
  coord_cartesian(ylim = c(0, 5), xlim = c(-3, 3)) +
  geom_text_repel(data = df %>% filter(feature %in% top10_features),
                  aes(label = genus),
                  max.overlaps = Inf, 
                  key_glyph = draw_key_point,
                  nudge_x = 2,       # Nudges top 10 to the right
                  direction = "y",   # Keeps vertical alignment
                  hjust = 0,                  
                  force = 12,           # Increase repulsion force
                  box.padding = 0.8,   # More spacing between labels
                  point.padding = 0.3) +       # Align text to the left of the callout line
  
  geom_text_repel(data = df %>% filter(feature %in% bottom10_features),
                  aes(label = genus),
                  max.overlaps = Inf,
                  key_glyph = draw_key_point,
                  nudge_x = -2,      # Nudges bottom 10 to the left
                  direction = "y",   # Keeps vertical alignment
                  hjust = 0, 
                  force = 12,           # Increase repulsion force
                  box.padding = 0.8,   # More spacing between labels
                  point.padding = 0.3)+  # More spacing from points# Align text to the right of the callout line
  theme(
    plot.margin = margin(20, 20, 20, 50),  # Adjust the margins (top, right, bottom, left)
    plot.title = element_text(size = rel(1.1),face="bold",hjust = 0.5)  # Make title font smaller and center it
  )

ggsave("volcano_genus.png", plot = volcano, width = 8, height = 6, dpi = 300)
