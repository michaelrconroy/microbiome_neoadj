rm(list=ls())
library(RColorBrewer)
library(ggplot2)
library(grid)
library(MASS)
library("pheatmap")
library(viridis)

# ------------------------------------------------------------------------
# set the working directory as the location of the script ---

setwd("S:/MolecularOncology/Users/Michael C/neoadjuvant microbiome/NGS results/microbiome_R_directory/")

# ------------------------------------------------------------------------
# create output directory
analysisdir = paste("S:/MolecularOncology/Users/Michael C/neoadjuvant microbiome/NGS results/microbiome_R_directory/analysis/heatmaps/genus_truncated", "/", sep="")


# ------------------------------------------------------------------------
inDir  = "indir/";
files  = c("jhmi-neoadj-mbiom-2024.genus.txt")

colscheme  = c("#ce4040",
               "#d68787",
               "#336699",
               "#84a3c1",
               "#F1C40F",
               "#27AE60",
               "#16A085",
               "#8E44AD",
               "#7F8C8D")

colscheme2 = c("#27AE60",
               "#ba87d0",
               "#F1C40F")

for (f in files){
  print(paste(f, "...", sep=""))
  tablefile = paste(inDir,f,sep="")
  A <- read.table(tablefile, sep="\t", header=FALSE, check.names=FALSE, as.is=TRUE)
  

  # Extract the first row for columns 29 to 259
  full_names <- A[1, 29:259]
  
  # Identify unassigned genera
  is_unass <- grepl("g_unassigned$", full_names)
  
  # Initialize output
  genus_names <- character(length(full_names))
  
  # Case 1: assigned genus → genus only
  genus_names[!is_unass] <-
    sub(".*g_", "", full_names[!is_unass])
  
  # Case 2: unassigned genus → f_<family>.g_unass
  genus_names[is_unass] <-
    paste0(
      "f_",
      sub(".*f_", "", sub("\\.g_.*", "", full_names[is_unass])),
      ".g_UA"
    )
  
  # Replace the first row
  A[1, 29:259] <- genus_names
  
  
  
  A = t(A)
  colnames(A)     = A[1,]
  rownames(A)     = A[,1]
  A               = A[2:nrow(A),2:ncol(A)]
  A               = A[,3:ncol(A)]  
  meta            = t(A[1:27,])
  meta            = data.frame(meta)
  meta            = meta[meta$any_PDL !=3, ]

  # ---
  # critical VARS ---
  meta$study                   = as.character(meta$study)
  meta$Record_ID               = as.character(meta$Record_ID)
  meta$Days_from_IO_start      = as.numeric(as.character(meta$Days_from_IO_start))
  meta$MPR                     = factor(meta$MPR, levels=c("0", "1"))
  meta$Response_6mo            = factor(meta$Response_6mo, levels=c("0", "1"))
  meta$percent_RVT             = factor(ifelse(as.numeric(gsub(">","",meta$percent_RVT)), "<=0.1", ">0.1"))
  meta$pretreatment_AJCC_stage = factor(meta$pretreatment_AJCC_stage)
  meta$Sex                     = factor(meta$Sex, levels=c("0", "1"))
  meta$smoking_status          = factor(meta$smoking_status)
  meta$age_enrollment          = as.numeric(as.character(meta$age_enrollment))
  meta$timept                  = ifelse(meta$Days_from_IO_start < 1, "pre-tx", "on-tx")
  meta$timept                  = factor(meta$timept , levels=c("pre-tx", "on-tx"))
  meta$any_PDL                 = factor(meta$any_PDL, levels=c("0","1","2"))

  # --------------------------------------------------------------------------
  # normalize values -----------
  A    = A[28:nrow(A),]
  A    = as.matrix(A)
  for (j in 1:nrow(A)){
    rownames(A)[j] = substr(rownames(A)[j],0,95)
  }
  all.vals = array(0, dim(A))
  for (i in 1:nrow(A)){
    for (j in 1:ncol(A)){
      all.vals[i,j] = as.numeric(as.character(A[i,j]))
    }
  }
  colnames(all.vals) = colnames(A)
  rownames(all.vals) = rownames(A)

  # --------------------------------------------------------------------------
  # color scheme -------------------------------------------------------------
  mycolors = c()

  mycolors$MPR["0"]   = "#ed5565"
  mycolors$MPR["1"]   = "#1ab394"

  mycolors$Response_6mo["0"]   = "#ec87c0"
  mycolors$Response_6mo["1"]   = "#48cfad"
  
  mycolors$any_PDL["0"]        = "gray92"
  mycolors$any_PDL["1"]        = "#4fc1e8"
  mycolors$any_PDL["2"]        = "#2980B9"

  mycolors$percent_RVT["<=0.1"]  = "#4fc1e8"
  mycolors$percent_RVT[">0.1"]   = "#ffce54"

  mycolors$Sex["0"]   = "#ff9d9d"
  mycolors$Sex["1"]   = "#83c5df"

  mycolors$study["J1414"]   = "#b17eb8"
  mycolors$study["J1772"]   = "#e6cae9"

  mycolors$timept["pre-tx"] = "#ffc99b"
  mycolors$timept["on-tx"]  = "#ff8f31"

  mycolors$smoking_status["0"]  = "#cccccc"
  mycolors$smoking_status["1"]  = "#ffc36a"
  mycolors$smoking_status["2"]  = "#da9551"

  # end color scheme ----------------------------------------------------------
  # --------------------------------------------------------------------------
  colscheme <- c("gray92", "#2C3E50", "#2980B9", "#8E44AD", "#C0392B", "#D35400", "#F39C12", "#F1C40F", "#27AE60", "#196f3e")
  CnormedCOLORS = colorRampPalette(colscheme[1:10])(20)
  CnormedCOLORS   = colorRampPalette(rev(brewer.pal(9,"YlGnBu")))(20)
  # end heatmap scheme ------------------------------------------------------------
  # visualize by subset of samples ----

  vals = log10((all.vals/100)+0.0001)
  vals = vals[,rownames(meta)]
  # select top ---
  vals = vals[order(rowMeans(vals), decreasing=TRUE),]
  selectN = min(nrow(vals),50)
  vals = vals[1:selectN,]
  meta = meta[order(meta$MPR, meta$Response_6mo, meta$any_PDL, meta$percent_RVT, meta$study, meta$timept, meta$smoking_status, meta$Sex),]
  meta = meta[,c("MPR", "Response_6mo", "percent_RVT", "any_PDL", "study", "timept", "smoking_status", "Sex")]
  vals = vals[,rownames(meta)]
  
  row_means <- rowMeans(vals)
  
  # Get top 10 genera by mean abundance
  top10 <- sort(row_means, decreasing = TRUE)[1:10]
  
  # Print nicely
  print(top10)
  

  pdf(file=paste(analysisdir,gsub("(jhmi-neoadj-mbiom-2024.|.txt)","",f),  ".clustered-euclidean.pdf", sep=""), width=9, height=7)
  pheatmap(vals, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=TRUE, annotation_col=meta, fontsize_row=4, fontsize = 7, annotation_colors=mycolors, scale="none", color=CnormedCOLORS, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean",
  border_color = "black")
  dev.off()

  pdf(file=paste(analysisdir,gsub("(jhmi-neoadj-mbiom-2024.|.txt)","",f),  ".ordered-euclidean.pdf", sep=""), width=9, height=7)
  pheatmap(vals, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=FALSE, annotation_col=meta, fontsize_row=4, fontsize = 7, annotation_colors=mycolors, scale="none", color=CnormedCOLORS, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean",
  border_color = "black")
  dev.off()

  pdf(file=paste(analysisdir,gsub("(jhmi-neoadj-mbiom-2024.|.txt)","",f),  ".clustered-manhattan.pdf", sep=""), width=9, height=7)
  pheatmap(vals, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=TRUE, annotation_col=meta, fontsize_row=4, fontsize = 7, annotation_colors=mycolors, scale="none",  color=CnormedCOLORS, clustering_distance_rows="manhattan", clustering_distance_cols="manhattan",
  border_color = "black")
  dev.off()

  pdf(file = paste0(analysisdir, gsub("(jhmi-neoadj-mbiom-2024.|.txt)", "", f), ".ordered-manhattan.pdf"),
      width = 9, height = 7)
  
  pheatmap(
    vals,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,   # must be TRUE to display labels
    cluster_cols = FALSE,
    annotation_col = meta,
    labels_col = meta$Record_ID,  # ← your Record_IDs
    angle_col = 90,               # rotate text vertically (bottom to top)
    fontsize_row = 4,
    fontsize = 7,
    annotation_colors = mycolors,
    scale = "none",
    color = CnormedCOLORS,
    clustering_distance_rows = "manhattan",
    clustering_distance_cols = "manhattan",
    border_color = "black"
  )
  
  dev.off()
}

# outputting a table of top 10 most abundant genera

# Assuming 'top10' is a named numeric vector like:
# Bacteroides = -0.9657898, Clostridium_XlVa = -1.1229043, etc.

# 1. Convert to a data frame
top10_df <- data.frame(
  Genus = names(top10),
  log10_value = as.numeric(top10)
)

# 2. Back-transform to percent abundance
top10_df$Mean_Abundance_Percent <- (10^(top10_df$log10_value) - 0.0001) * 100

# 3. Round for presentation
top10_df$Mean_Abundance_Percent <- round(top10_df$Mean_Abundance_Percent, 2)
top10_df$log10_value <- round(top10_df$log10_value, 3)

# 4. Order nicely by decreasing abundance
top10_df <- top10_df[order(-top10_df$Mean_Abundance_Percent), ]

# 5. Save as Excel file
library(openxlsx)
write.xlsx(top10_df, file = "top10_genera_mean_abundance.xlsx", row.names = FALSE)

