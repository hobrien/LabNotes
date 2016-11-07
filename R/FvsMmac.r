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

workDir <- "~/BTSync/FetalRNAseq/Counts/MvsF_14_20_noA_excl_15641_18432_Cooks.75"      # working directory for the R session

projectName <- "MvsF_14_20_noA_excl_15641_18432_Cooks.75"                         # name of the project
author <- "Heath O'Brien"                                # author of the statistical analysis/report

rawDir <- "~/BTSync/FetalRNAseq/Counts/raw"                                      # path to the directory containing raw counts files
featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")# NULL if no feature to remove

varInt <- "Sex"                                    # factor of interest
condRef <- "Female"                                      # reference biological condition
batch <- c("PCW", "Centre", "RIN")                # blocking factor: NULL (default) or "batch" for example
interact <- NULL #c("PCW")
RIN_cutoff <- 0
PCW_cutoff <- c(14, 20)
fitType <- "parametric"                              # mean-variance relationship: "parametric" (default) or "local"
cooksCutoff <- 100000                             # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.05                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"
testMethod <- 'Wald'
typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors

colors <- c("dodgerblue","firebrick1",               # vector of colors of each biological condition on the plots
            "MediumVioletRed","SpringGreen")
exclude <- c('15641', '18432', '16428')
################################################################################
###                             running script                               ###
################################################################################
dir.create(workDir)
setwd(workDir)
library(devtools)
load_all(pkg = "~/BTSync/Code/R/SARTools")
library(dplyr)
library(tidyr)
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
target <- filter(target, ! grepl('A', label))
#target <- mutate(target, PCW = factor(floor(PCW)))
# loading counts
# this doesn't use the batch info
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# description plots
# this doesn't use the batch info
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
if (testMethod=='Wald' ) {
  out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch, interact=interact,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)
} else if (testMethod=='LRT' ) {
  out.DESeq2 <- run.DESeq2.LRT(counts=counts, target=target, varInt=varInt, batch=batch, interact=interact,
                           locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                           cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)
  } else {
  stop("testMethod nor recognised")
}
mcols(out.DESeq2$dds)$maxCooks <- apply(assays(out.DESeq2$dds)[["cooks"]], 1, max)
out.DESeq2$results$Male_vs_Female$pvalue[mcols(out.DESeq2$dds)$maxCooks > cooksCutoff] <- NA
out.DESeq2$results$Male_vs_Female$padj <- p.adjust(out.DESeq2$results$Male_vs_Female$pvalue, method="BH")

# PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)

# Filter 

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
MalevsFemale.up <- read.delim("tables/MalevsFemale.up.txt", check.names=FALSE)
MalevsFemale.up <- filter(MalevsFemale.up, padj < 0.05) 
MalevsFemale.up <- select(MalevsFemale.up, Id, Female, Male, FoldChange, log2FoldChange, pvalue, padj)
MalevsFemale.up <- bind_cols(GetGeneIDs(MalevsFemale.up$Id), MalevsFemale.up)
write.table(MalevsFemale.up, file="tables/MaleUp.txt", sep="\t", quote=FALSE, row.names=FALSE)

MalevsFemale.down <- read.delim("tables/MalevsFemale.down.txt", check.names=FALSE)
MalevsFemale.down <- filter(MalevsFemale.down, padj < 0.05) %>% select(Id, Female, Male, FoldChange, log2FoldChange, pvalue, padj)
MalevsFemale.down <- bind_cols(GetGeneIDs(MalevsFemale.down$Id), MalevsFemale.down)
write.table(MalevsFemale.down, file="tables/FemaleUp.txt", sep="\t", quote=FALSE, row.names=FALSE)

# This produces slightly different numbers than the independent filtering used by DESeq2
# (eg; filtering 40298 vs. 40342). I assume this has something to do with rounding error
# when means are calculated?
MalevsFemale.complete <- read.delim("tables/MalevsFemale.complete.txt")
MalevsFemale.complete$CountMean <- select(MalevsFemale.complete, starts_with('norm')) %>% rowMeans()
MalevsFemale.complete <- filter(MalevsFemale.complete, CountMean >= as.numeric(tabIndepFiltering(out.DESeq2$results)[2]))

MalevsFemale.complete <- bind_cols(GetGeneIDs(MalevsFemale.complete$Id), MalevsFemale.complete)
MalevsFemale.complete <-  separate(MalevsFemale.complete, Id, c("Id"), sep='[.]', extra='drop')
select(MalevsFemale.complete, Id, GeneID,Female,	Male,	FoldChange,	log2FoldChange,	pvalue,	padj) %>% 
  write.table(file="tables/Background2.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Make gtc file for GSEA
write('#1.2', file = "tables/MvsF.gct")
MalevsFemale.gct <- select(MalevsFemale.complete, NAME=GeneID, DESCRIPTION=Id, starts_with('norm'))
colnames(MalevsFemale.gct) <- gsub("norm.", "", colnames(MalevsFemale.gct))
write(c(nrow(MalevsFemale.gct), ncol(MalevsFemale.gct)-2), file = "tables/MvsF.gct",
      #ncolumns = if(is.character(x)) 1 else 5,
      append = TRUE, sep = "\t")
write.table(MalevsFemale.gct, file="tables/MvsF.gct", sep="\t", quote=FALSE, row.names=FALSE, append=TRUE)

# Make cls file for GSEA
write(paste(nrow(target), 2, 1), file = "tables/MvsF.cls")
write("# Female Male", file = "tables/MvsF.cls", append=TRUE)
write(paste(target$Sex, collapse= " "), file = "tables/MvsF.cls", append=TRUE)

# Include histograms of PCW and RIN
ggplot(target, aes(x=PCW)) +
  geom_histogram(binwidth=1/7) +
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
DE.cooks <- as.data.frame(assays(out.DESeq2$dds)[["cooks"]][as.character(rbind(MalevsFemale.up, MalevsFemale.down)$Id),])
gather(DE.cooks, 'Sample', 'Cooks')  %>% mutate(Cooks = ifelse(Cooks > 100, 100, Cooks)) %>%
  ggplot(aes(x=Cooks))+geom_histogram(bins=100)+ facet_wrap(~Sample) + coord_cartesian(ylim=c(0, 10))
ggsave("figures/CooksHist.png")
