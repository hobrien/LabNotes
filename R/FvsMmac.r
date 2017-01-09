################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Hugo Varet
### May 9th, 2016
### designed to be executed with SARTools 1.3.0
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
rm(list=ls())                                        # remove all the objects from the R session

projectName <- "MvsF_12_14_Cooks.75_excl_16491_FDR.1_new"                         # name of the project

workDir <- paste("~/BTSync/FetalRNAseq/Counts", projectName, sep='/')      # working directory for the R session

author <- "Heath O'Brien"                                # author of the statistical analysis/report

rawDir <- "~/BTSync/FetalRNAseq/Counts/raw"                                      # path to the directory containing raw counts files
featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")# NULL if no feature to remove

varInt <- "Sex"                                    # factor of interest
condRef <- "Female"                                      # reference biological condition
batch <- c("PCW", "Centre", "RIN")                # blocking factor: NULL (default) or "batch" for example
interact <- c()
RIN_cutoff <- 0
PCW_cutoff <- c(12, 14)
fitType <- "parametric"                              # mean-variance relationship: "parametric" (default) or "local"
#if numeric, features with maxCooks values above this number are removed 
cooksCutoff <-  FALSE                          # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.1                                    # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"
testMethod <- 'Wald'
typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors

BrainBank <- 'HDBR' # 'All' #
exclude <- c()
exclude <- c('16491')
#exclude <- c('15641', '18432')#, '16491')
#exclude <- c("15641", "16548", "17160", "17923", "18294", "18983", "17921", "17486", "16024", "16115", "16810", "16826", "17048", "17053", "17071", "17333", "18432", "18666", "17264")
colors <- c("dodgerblue","firebrick1",               # vector of colors of each biological condition on the plots
            "MediumVioletRed","SpringGreen")
excludedFeaturesFile = NA
#excludedFeaturesFile = "~/BTSync/FetalRNAseq/LabNotes/SexChrGenes.txt" 
################################################################################
###                             running script                               ###
################################################################################
dir.create(workDir)
setwd(workDir)
library(devtools)
load_all(pkg = "~/BTSync/Code/R/SARTools")
library(dplyr)
library(tidyr)
library(readr)
source("~/BTSync/FetalRNAseq/LabNotes/R/FormatGGplot.R")

# checking parameters
# this checks that there is a single batch factor
#checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
#                       rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
#                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
#                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
#                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)

# loading target file
# this checks the batch so it will need to be modified
# once this is loaded, I can filter it to remove samples with RIN < 5
#target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)
targetFile <- "~/BTSync/FetalRNAseq/LabNotes/MvsFmac.txt"
target <- read.delim(targetFile)                        # path to the design/target file

sample_info <- read.delim("~/BTSync/FetalRNAseq/LabNotes/sample_info.txt")
target <- left_join(target, select(sample_info, BrainBankID, Sex, PCW, RIN), by = c("label" = "BrainBankID"))

sample_progress <- read.delim("~/BTSync/FetalRNAseq/LabNotes/SampleProgress.txt")
target <- left_join(target, select(sample_progress, sample, Centre), by = c("label" = "sample"))
target <- arrange(target, Sex)
if (!is.null(RIN_cutoff)) {
  target <- filter(target, RIN >= RIN_cutoff)
}
if (!is.null(PCW_cutoff)) {
  target <- filter(target, PCW >= PCW_cutoff[1] & PCW < PCW_cutoff[2])
}
if (length(exclude) > 0) {
  target <- filter(target, !label %in% exclude)
}
if (BrainBank == 'HDBR') {
    target <- filter(target, ! grepl('A', label))
}
target <- mutate(target, PCW = floor(PCW))
# loading counts
# this doesn't use the batch info
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)
if (! is.na(excludedFeaturesFile)) {
  excludedFeatures <- read_csv(excludedFeaturesFile, col_names = FALSE)
  counts <- counts[!rownames(counts) %in% excludedFeatures$X1, ]
}
# description plots
# this doesn't use the batch info
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
if (testMethod=='Wald' ) {
  out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch, interact=interact,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)
  mcols(out.DESeq2$dds)$maxCooks <- apply(assays(out.DESeq2$dds)[["cooks"]], 1, max)
  if (is.numeric(cooksCutoff)) {
    out.DESeq2$results$Male_vs_Female$pvalue[mcols(out.DESeq2$dds)$maxCooks > cooksCutoff] <- NA
    out.DESeq2$results$Male_vs_Female$padj <- p.adjust(out.DESeq2$results$Male_vs_Female$pvalue, method="BH")
  }  
} else if (testMethod=='LRT' ) {
  out.DESeq2 <- run.DESeq2.LRT(counts=counts, target=target, varInt=varInt, batch=batch, interact=interact,
                           locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                           cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)
} else {
stop("testMethod not recognised")
}

# PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)

#Save VST counts
vst <- as.data.frame(assay(varianceStabilizingTransformation(out.DESeq2$dds)))
vst$Id <- row.names(vst)
vst <- vst[,c(ncol(vst), 1:ncol(vst)-1)]
write.table(vst, file="tables/VST.txt", sep="\t", quote=FALSE, row.names = FALSE)

#filter out genes with counts < 10 in more than 10% of samples
filterSet <- rowSums(ifelse(counts(out.DESeq2$dds) > 9, 1, 0)) > (ncol(vst)-1)*.9
write.table(vst[filterSet, ], file="tables/VST_filtered.txt", sep="\t", quote=FALSE, row.names = FALSE)

filterSet <- rowSums(ifelse(counts(out.DESeq2$dds) > 9, 1, 0)) > (ncol(vst)-1)*.5
write.table(vst[filterSet, ], file="tables/VST_filtered2.txt", sep="\t", quote=FALSE, row.names = FALSE)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)





my_db <- src_mysql("FetalRNAseq", host="localhost", user="root")
GetGeneIDs <- function(Ids) {
  geneIDs <- tibble(Chr=character(), GeneID=character())
  for (Id in Ids) {
  geneIDs <- collect(tbl(my_db, sql(paste0("SELECT DISTINCT GencodeGTF.seqid AS 'Chr', 
                                   GencodeFeatures.value AS 'GeneID' FROM 
                                   GencodeFeatures, GencodeGTF WHERE 
                                   GencodeGTF.id = GencodeFeatures.id AND 
                                   GencodeFeatures.feature = 'gene_name' AND 
                                   GencodeFeatures.id IN (
                                   SELECT GencodeGTF.id FROM GencodeGTF, GencodeFeatures WHERE 
                                   GencodeGTF.id = GencodeFeatures.id AND GencodeFeatures.Value =",
                                paste0("'", Id, "'"),
                                ")"
                                )))) %>%
    bind_rows(geneIDs, .)
  }
  geneIDs
}
DEgenes=NA
if (testMethod == 'Wald') {
  MalevsFemale.up <- read.delim("tables/MalevsFemale.up.txt", check.names=FALSE)
  MalevsFemale.up <- select(MalevsFemale.up, Id, Female, Male, FoldChange, log2FoldChange, pvalue, padj)
  MalevsFemale.up <- bind_cols(GetGeneIDs(MalevsFemale.up$Id), MalevsFemale.up)
  write.table(MalevsFemale.up, file="tables/MaleUp.txt", sep="\t", quote=FALSE, row.names=FALSE)
  MalevsFemale.down <- read.delim("tables/MalevsFemale.down.txt", check.names=FALSE)
  MalevsFemale.down <- select(MalevsFemale.down, Id, Female, Male, FoldChange, log2FoldChange, pvalue, padj)
  MalevsFemale.down <- bind_cols(GetGeneIDs(MalevsFemale.down$Id), MalevsFemale.down)
  write.table(MalevsFemale.down, file="tables/FemaleUp.txt", sep="\t", quote=FALSE, row.names=FALSE)
  DEgenes <- nrow(MalevsFemale.down) + nrow(MalevsFemale.up)
  DESeq.complete <- read.delim("tables/MalevsFemale.complete.txt")
} else if (testMethod=='LRT' ) {
  lrt.up <- read.delim(paste0("tables/drop", varInt, ".up.txt"), check.names=FALSE)
  lrt.up <- select(lrt.up, Id, baseMean, FoldChange, log2FoldChange, pvalue, padj)
  lrt.up <- bind_cols(GetGeneIDs(lrt.up$Id), lrt.up)
  write.table(lrt.up, file=paste0("tables/drop", varInt, ".up.gene_name.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  lrt.down <- read.delim(paste0("tables/drop", varInt, ".down.txt"), check.names=FALSE)
  lrt.down <- select(lrt.down, Id, baseMean, FoldChange, log2FoldChange, pvalue, padj)
  lrt.down <- bind_cols(GetGeneIDs(lrt.down$Id), lrt.down)
  write.table(lrt.down, paste0("tables/drop", varInt, ".down.gene_name.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  DEgenes <- nrow(lrt.up) + nrow(lrt.down)
  DESeq.complete <- read.delim(paste0("tables/drop", varInt, ".down.txt"))
} else {
  stop("testMethod not recognised")
}

# This produces slightly different numbers than the independent filtering used by DESeq2
# (eg; filtering 40298 vs. 40342). I assume this has something to do with rounding error
# when means are calculated?
DESeq.complete$CountMean <- select(DESeq.complete, starts_with('norm')) %>% rowMeans()
DESeq.complete <- filter(DESeq.complete, CountMean >= as.numeric(tabIndepFiltering(out.DESeq2$results)[2]))

#MalevsFemale.complete <- bind_cols(GetGeneIDs(MalevsFemale.complete$Id), MalevsFemale.complete)
DESeq.complete <-  separate(DESeq.complete, Id, c("Id"), sep='[.]', extra='drop')
select(DESeq.complete, Id, Female,	Male,	FoldChange,	log2FoldChange,	pvalue,	padj) %>% 
  write.table(file="tables/Background.txt", sep="\t", quote=FALSE, row.names=FALSE)

#write summary of analysis to file
summary <- data.frame(BrainBank=c(BrainBank), 
                      AgeRange=c(paste(PCW_cutoff, collapse='-')), 
                      RIN=c(ifelse(RIN_cutoff==0, "All", RIN_cutoff)),
                      FDR=c(alpha), 
                      Excluded=c(ifelse(length(exclude)==0,
                                        "None",
                                        ifelse(length(exclude)>3,
                                               length(exclude),
                                               paste(exclude, collapse="/")
                                              )
                                       )
                                ),
                      n=c(ncol(counts)),
                      test=c(testMethod),
                      model=c(ifelse(length(interact)==0, '+', '*')), 
                      CooksCutoff=c(ifelse(testMethod=='Wald', cooksCutoff, "None")),
                      DEgenes=c(DEgenes),
                      res=c(workDir)                     
                      )                                              
write.table(summary, 
            file="~/BTSync/FetalRNAseq/Counts/Summary.txt", 
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE, 
            col.names=FALSE,
            append=TRUE
            )


# Include histograms of PCW and RIN
ggplot(target, aes(x=PCW)) +
  geom_bar() +
  facet_grid(Sex ~ .) +
  tufte_theme() +
  scale_y_continuous(breaks=seq(0,100, 2))
ggsave("figures/PCW_hist.png")

ggplot(target, aes(x=RIN)) +
  geom_histogram(binwidth=1) +
  facet_grid(Sex ~ .) +
  scale_y_continuous(breaks=seq(0,100,2)) +
  scale_x_continuous(breaks=seq(0,12,1)) +
  tufte_theme()
ggsave("figures/RIN_hist.png")

# Plot histogram of Cooks distances for DE genes
# last bin should be read as '100+')
# counts of 10 should be read as '10+'
#this only runs when cooksCutoff is set to NA becasue it doesn't make sense when high maxCooks features are excluded
if (testMethod == 'Wald' & ! is.numeric(cooksCutoff)) {
  DE.cooks <- as.data.frame(assays(out.DESeq2$dds)[["cooks"]][as.character(rbind(MalevsFemale.up, MalevsFemale.down)$Id),])
  gather(DE.cooks, 'Sample', 'Cooks')  %>% mutate(Cooks = ifelse(Cooks > 100, 100, Cooks)) %>%
    ggplot(aes(x=Cooks)) +
      geom_histogram(bins=100) + 
      facet_wrap(~Sample) + 
      coord_cartesian(ylim=c(0, 10)) +
      tufte_theme() +
      theme(axis.text.x=element_text(size=6)) +
      theme(axis.text.y=element_text(size=6)) +
      theme(strip.text=element_text(size=4))
  
  ggsave("figures/CooksHist.png")
}
