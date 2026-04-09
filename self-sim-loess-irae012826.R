rm(list=ls())
library(RColorBrewer)
library(ggplot2)
library(grid)
library(MASS)
library("pheatmap")
library(vegan)
library(ape)
require(GGally)
require(gtools)
library(ggpubr)
library(reshape2)
library(readxl)
# set the working directory to the exact code base location ------------------
setwd("S:/MolecularOncology/Users/Michael C/neoadjuvant microbiome/NGS results/microbiome_R_directory")
# set the place for the analysis outputs ------------------
analysisdir = "analysis/irAE/LOESS/"
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)
set.seed(5432)

# load data from file in A01.1-add-metadata.pl
alphaDf <- read_excel("jhmi-neoadj-mbiom-2024.alpha-diversity.xlsx", col_names = TRUE)

# load data from file
meta <- read.table("shared-species-prct-results.txt", header=TRUE, sep="\t")

# limit first and second column to matching ids in alphaDf$SampleID
meta          = meta[meta$SimpleID.1 %in% alphaDf$SampleID, ]
meta          = meta[meta$SimpleID.2 %in% alphaDf$SampleID, ]

# format variables
meta$timept.1 = as.character(meta$timept.1)
meta$timept.2 = as.character(meta$timept.2)

meta$timept   = as.numeric(as.character(meta$timept.1))
meta$MPR      = as.numeric(as.character(meta$MPR.1))


meta$comp.type = as.character(meta$comp.type)

# Self vs non-self check
meta$compcat       = "self"

# drop non-shared reads similarity measures
meta = meta[meta$dist.type == "shared.species.prct",]

# timepoint limitations
# meta            = meta[meta$timept >= -20,]
meta            = meta[meta$timept <= 365,]

# Order by timepoint variable
meta            = meta[order(meta$timept),]

# Exclude perfect zero distance matches (if present)
# meta            = meta[!(meta$dist == 0), ]

# ---
# Convert negative values days to 0
meta[meta$timept < 0, "timept"] <- 0
# ----
# add meta vars
metaA          = read.table("meta.tsv", sep="\t", header=TRUE, check.names=FALSE)
metaA          = unique(metaA[,c("Record_ID", "pretreatment_AJCC_stage", "Sex", "Tumor", "age_enrollment", "smoking_status",
                                 "tx_type", "MPR", "percent_RVT", "Response_6mo", "any_irae", "high_grade_irae")])
# ----
meta$Record_ID = meta$Record_ID.1
meta           = dplyr::left_join(meta, metaA, by="Record_ID")
meta$MPR       = meta$MPR.y
# ----
meta$Sex       = ifelse(meta$Sex=="1",   "male", "female")
meta$Tumor     = ifelse(meta$Tumor=="1", "AD", "SQ")
# meta$MPR       = ifelse(meta$MPR=="1", "mpr", "nr")
meta$any_irae <- as.factor(meta$any_irae)

# nonparametric difference test
npdiff <- function(a,b){
  a = a[!is.na(a)]
  b = b[!is.na(b)]
  
  perms    = 1000
  diff0    = abs(mean(a)-mean(b))
  na       = length(a)
  nb       = length(b)
  
  if (length(a) == 0){
    return(1)
  }
  if (length(b) == 0){
    return(1)
  }
  combined = c(a,b)
  simdiffs = c()
  for (j in 1:perms){
    ncomb    = sample(combined)
    newa     = ncomb[1:na]
    newb     = ncomb[(na+1):length(ncomb)]
    simdiffs = c(simdiffs, abs(mean(newa)-mean(newb)))
  }
  simdiffs   = simdiffs[order(simdiffs, decreasing=TRUE)]
  mymaxdex = 0
  for (j in 1:perms){
    if (diff0 < simdiffs[j]){
      # stop and record
      mymaxdex = j
    }
  }
  if (mymaxdex==0){
    mymaxdex = 1
  }
  inferredP = mymaxdex/perms
  return(inferredP)
}


# ---
# Check to confirm 1 d0 per patient
# table(meta[meta$timept == 0, "Record_ID.1"])
if (max(table(meta[meta$timept == 0, "Record_ID.1"])) > 1){
  stop("Error too many baseline samples by Record_ID")
}
if (max(table(meta[meta$timept == 0, "Record_ID.2"])) > 1){
  stop("Error too many baseline samples by Record_ID")
}

# colors ---
# colors ---
mycolors = list()
mycolors$MPR["0"]  = "#f94144"
mycolors$MPR["1"]  = "#277da1"


mycolors$any_irae["0"]          = "#ff7f0e"
mycolors$any_irae["1"]          = "#2ca02c"

meta$MPR           = factor(meta$MPR)

dists = c("shared.species.prct")
for (mydist in dists[1]){
  p1 <- ggplot(meta, aes(x=timept,y=dist)) +
    geom_path(aes(group=Record_ID.1, color=any_irae), alpha=0.2) +
    geom_point(aes(color=any_irae), pch=19, size=2.5, alpha=0.4) +
    scale_fill_manual(values=mycolors$any_irae)  +
    scale_color_manual(values=mycolors$any_irae) +
    xlab("Days from IO Baseline") +
    ylab("Distance Measure") +
    theme_bw() +
    theme(axis.text.x  = element_text(size=9, colour="black"),
          axis.text.y  = element_text(size=9, colour="black"),
          axis.title.x = element_text(size=9, colour="black"),
          axis.title.y = element_text(size=9, colour="black"),
          plot.title   = element_text(size=10, colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio=0.75,
          legend.title = element_blank()) +
    ggtitle(paste0(mydist, " all samples self by "))
  ggsave(paste(analysisdir, mydist, "-00.selected.pdf", sep=""), p1, width=8, height=5)
  write.csv(meta[meta$dist.type==mydist,], file=paste(analysisdir, mydist,"-00.selected.csv", sep=""), row.names=FALSE)
  


  # R vs NR loess bootstrap
  rDat  = meta[meta$dist.type==mydist & meta$any_irae == "1",]
  nrDat = meta[meta$dist.type==mydist & meta$any_irae == "0",]
  
  # ---
  # get bootstrapped loess levels at 30, 60, 90, and 180
  bootStats  = c()
  bootProp   = 0.40
  for (subi in 1:10){
    subi_r_dat = rDat[sample(1:nrow(rDat))[1:floor(nrow(rDat)*bootProp)],]
    while(min(subi_r_dat$timept) > 30 | max(subi_r_dat$timept) < 90){
      subi_r_dat = rDat[sample(1:nrow(rDat))[1:floor(nrow(rDat)*bootProp)],]
    }
    
    rDatMi   <- loess(dist~timept, subi_r_dat, span=10)
    bootPredr = predict(rDatMi, newdata = data.frame(timept=c(30,60,90,180)))
    
    
    subi_nr_dat = nrDat[sample(1:nrow(nrDat))[1:floor(nrow(nrDat)*bootProp)],]
    while(min(subi_nr_dat$timept) > 30 | max(subi_nr_dat$timept) < 90){
      subi_nr_dat = nrDat[sample(1:nrow(nrDat))[1:floor(nrow(nrDat)*bootProp)],]
    }
    
    nrDatMi   <- loess(dist~timept, subi_nr_dat, span=10)
    bootPrednr = predict(nrDatMi, newdata = data.frame(timept=c(30,60,90,180)))
    
    bootStats = rbind(bootStats,c(subi,bootPredr,bootPrednr))
  }
  bootStats = data.frame(bootStats)
  colnames(bootStats) = c("subi", "r.d30",  "r.d60",  "r.d90",  "r.d180",
                          "nr.d30", "nr.d60", "nr.d90", "nr.d180")
  
  write.csv(bootStats, file=paste(analysisdir, mydist, ".1.03l-d180-time-series.csv", sep=""), row.names=FALSE)
  
  thismnR30  = sprintf("%3.2f", mean(bootStats$r.d30, na.rm=TRUE))
  thismnnR30 = sprintf("%3.2f", mean(bootStats$nr.d30, na.rm=TRUE))
  thisMWP30  = sprintf("%3.3f", npdiff(bootStats$r.d30, bootStats$nr.d30))
  
  thismnR60  = sprintf("%3.2f", mean(bootStats$r.d60, na.rm=TRUE))
  thismnnR60 = sprintf("%3.2f", mean(bootStats$nr.d60, na.rm=TRUE))
  thisMWP60  = sprintf("%3.3f", npdiff(bootStats$r.d60, bootStats$nr.d60))
  
  thismnR90  = sprintf("%3.2f", mean(bootStats$r.d90, na.rm=TRUE))
  thismnnR90 = sprintf("%3.2f", mean(bootStats$nr.d90, na.rm=TRUE))
  thisMWP90  = sprintf("%3.3f", npdiff(bootStats$r.d90, bootStats$nr.d90))
  
  thismnR180  = sprintf("%3.2f", mean(bootStats$r.d180, na.rm=TRUE))
  thismnnR180 = sprintf("%3.2f", mean(bootStats$nr.d180, na.rm=TRUE))
  thisMWP180  = sprintf("%3.3f", npdiff(bootStats$r.d180, bootStats$nr.d180))
  
  myTitle = paste0("Bootstrap analysis results:\n",
                   "Day 30:  R mean: ",  thismnR30, "; NR mean: ", thismnnR30, "; Boot-Pval: ", thisMWP30, "\n",
                   "Day 60:  R mean: ",  thismnR60, "; NR mean: ", thismnnR60, "; Boot-Pval: ", thisMWP60, "\n",
                   "Day 90:  R mean: ",  thismnR90, "; NR mean: ", thismnnR90, "; Boot-Pval: ", thisMWP90, "\n",
                   "Day 180: R mean: ",  thismnR180, "; NR mean: ", thismnnR180, "; Boot-Pval: ", thisMWP180)
  
  p1 <- ggplot(meta[meta$dist.type==mydist,], aes(x=timept,y=dist)) +
    geom_vline(xintercept=c(30,60,90), color="black", linetype="dashed", alpha=0.95, size=0.2) +
    geom_smooth(data=rDat,  method = 'loess', color="#2ca02c", fill="#2ca02c",  span=10, alpha=0.25) +
    geom_smooth(data=nrDat, method = 'loess', color="#ff7f0e", fill="#ff7f0e", span=10, alpha=0.25) +
    xlab("Days from IO Baseline") +
    ylab("Self-similarity") +
    theme_bw() +
    theme(axis.text.x  = element_text(size=11, colour="black"),
          axis.text.y  = element_text(size=11, colour="black"),
          axis.title.x = element_text(size=11, colour="black"),
          axis.title.y = element_text(size=11, colour="black"),
          plot.title   = element_text(size=10, colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio=0.75,
          legend.title = element_blank()) +
    xlim(c(0,180)) +
    ggtitle(myTitle) +
    coord_cartesian(ylim=c(75,100))
  ggsave(paste(analysisdir, mydist, ".1.03l-d180-time-series.pdf", sep=""), p1, width=6, height=6)
  
  bootStatsMelt = melt(bootStats, id.vars=c("subi"))
  bootStatsMelt = bootStatsMelt[!grepl("180", bootStatsMelt$variable),]
  bootStatsMelt$any_irae = toupper(gsub("\\.d.+$", "", bootStatsMelt$variable))
  bootStatsMelt$timept        = paste0("Day ", gsub(".+\\.d", "", bootStatsMelt$variable))
  bootStatsMelt$timept        = factor(bootStatsMelt$timept, levels=c("Day 30", "Day 60", "Day 90"))
  
  p1 <- ggplot(bootStatsMelt, aes(x=any_irae,y=value)) +
    geom_boxplot(aes(color=any_irae), fill=NA, outlier.size=0, outlier.shape=NA, coef=1e100, alpha=0.6) +
    geom_point(aes(fill=any_irae), pch=21, color="black", alpha=0.6, size=2.5) +
    scale_fill_manual(values=mycolors$any_irae)  +
    scale_color_manual(values=mycolors$any_irae) +
    xlab(NULL) +
    ylab("Bootstrap Model Fit") +
    theme_bw() +
    theme(axis.text.x  = element_text(size=11, colour="black"),
          axis.text.y  = element_text(size=11, colour="black"),
          axis.title.x = element_text(size=11, colour="black"),
          axis.title.y = element_text(size=11, colour="black"),
          plot.title   = element_text(size=10, colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank()) +
    ggtitle(myTitle) +
    facet_wrap(~timept, ncol=3)
  ggsave(paste(analysisdir, mydist, ".1.03b-d180-time-series.pdf", sep=""), p1, width=6, height=6)
  


}
