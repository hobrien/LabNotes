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

workDir <- "/c8000xd3/rnaseq-heath/Counts/MvsF_all"      # working directory for the R session

projectName <- "MvsF"                         # name of the project
author <- "Heath O'Brien"                                # author of the statistical analysis/report

rawDir <- "/"                                      # path to the directory containing raw counts files
featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")# NULL if no feature to remove

varInt <- "treatment"                                    # factor of interest
condRef <- "F"                                      # reference biological condition
batch <- c("centre", "PCW", "RIN")                  # blocking factor: NULL (default) or "batch" for example
RIN_cutoff <- 0
fitType <- "parametric"                              # mean-variance relationship: "parametric" (default) or "local"
cooksCutoff <- TRUE                                  # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.1                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors

colors <- c("dodgerblue","firebrick1",               # vector of colors of each biological condition on the plots
            "MediumVioletRed","SpringGreen")

################################################################################
###                             running script                               ###
################################################################################
setwd(workDir)
library(devtools)
load_all(pkg = "~/src/SARTools")
#load_all(pkg = "~/BTSync/Code/R/SARTools")
library(dplyr)

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

target <- read.delim("~/LabNotes/MvsF.txt")                        # path to the design/target file
#target <- read.delim("~/BTSync/FetalRNAseq/LabNotes/MvsF.txt")                        # path to the design/target file

sample_info <- read.delim("~/LabNotes/sample_info.txt")
#sample_info <- read.delim("~/BTSync/FetalRNAseq/LabNotes/sample_info.txt")
target <- left_join(target, select(sample_info, BrainBankID, Sex, PCW, RIN), by = c("label" = "BrainBankID"))

sample_progress <- read.delim("~/LabNotes/SampleProgress.txt")
#sample_progress <- read.delim("~/BTSync/FetalRNAseq/LabNotes/SampleProgress.txt")
target <- left_join(target, select(sample_progress, sample, Centre), by = c("label" = "sample")) %>% View()

if (!is.null(RIN_cutoff)) {
  target <- filter(target, RIN >= RIN_cutoff)
}
# loading counts
# this doesn't use the batch info
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# description plots
# this doesn't use the batch info
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

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

