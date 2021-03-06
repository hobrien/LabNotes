################################################################################
### R script to compare several conditions with the SARTools and edgeR packages
### Hugo Varet
### May 9th, 2016
### designed to be executed with SARTools 1.3.0
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
rm(list=ls())                                        # remove all the objects from the R session
library("optparse")

option_list <- list(
  make_option(c("-m", "--min"), type="integer", default=NULL, 
              help="minimum age (PCW)", metavar="minimum age"),
  make_option(c("-x", "--max"), type="integer", default=NULL, 
              help="maximum age (PCW)", metavar="maximum age"),
  make_option(c("-r", "--rin"), type="numeric", default=NA, 
              help="minimum RIN", metavar="minimum RIN"),
  make_option(c("-p", "--pvalue"), type="numeric", default=0.1, 
              help="corrected pvalue cutoff", metavar="pvalue"),
  make_option(c("-b", "--brainbank"), type="character", default='HDBR', 
              help="Which brainbank to include ('HDBR', 'All', HDBRexpression)", metavar="brainbank"),
  make_option(c("-e", "--exclude"), type="character", default='15641,18432,16491', 
              help="Samples to exclude (comma separated list, no spaces)", metavar="excluded"),
  make_option(c("-s", "--sex_chromosomes"), action='store_true', type="logical", default=FALSE, 
              help="Exclude sex chromosomes", metavar="sex_chromosomes"),
  make_option(c("-a", "--age"), action='store_true', type="logical", default=FALSE, 
              help="Include age as a cofactor", metavar="age")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
#opt$min=14
#opt$max=15
#opt$brainbank='HDBRexpression'
if (is.null(opt$min) | is.null(opt$max)){
  print_help(opt_parser)
  stop("Min and Max ages must be specified", call.=FALSE)
}
if (opt$brainbank != 'HDBR' & opt$brainbank != 'All' & opt$brainbank != 'HDBRexpression'){
  print_help(opt_parser)
  stop("BrainBank options are 'HDBR', 'HDBRexpression' and 'All'", call.=FALSE)
}
if ( opt$max-opt$min < 1 ) {
  print_help(opt_parser)
  stop("Min age must be less than max", call.=FALSE)
}

PCW_cutoff <- c(opt$min, opt$max)
#PCW_cutoff <-c(12,20)
RIN_cutoff <- opt$rin
alpha <- opt$pvalue                                      # threshold of statistical significance
BrainBank <- opt$brainbank
pcw <- opt$age
exclude_sex <- opt$sex_chromosomes

exclude <- strsplit(opt$exclude, ',')[[1]]
if (BrainBank == 'HDBRexpression') {
  exclude <- exclude[exclude != "15641" & exclude != "18432" & exclude != "16491"]
} else {
  if (! findInterval(15, PCW_cutoff) == 1) {
    exclude <- exclude[exclude != "15641"]
  }
  if (! findInterval(17, PCW_cutoff) == 1) {
    exclude <- exclude[exclude != "18432"]
  }
  if (! findInterval(13, PCW_cutoff) == 1) {
    exclude <- exclude[exclude != "16491"]
  }
}  


projectName <- paste("MvsF", 
                     ifelse(PCW_cutoff[2]-PCW_cutoff[1] > 1, 
                            paste(PCW_cutoff, 
                                  collapse ='_'
                            ),
                            PCW_cutoff[1]
                     ),
                     ifelse(length(exclude > 0),
                            paste(c(BrainBank, 'excl', exclude), collapse='_', sep='_'),
                            BrainBank
                     ),
                     ifelse(pcw, "PCW_FDR", "FDR"),
                     alpha,
                     ifelse(exclude_sex, "noSexChr_edgeR", "edgeR"), 
                     sep='_'
)                         # name of the project


if (exclude_sex) {
  excludedFeaturesFile = "~/BTSync/FetalRNAseq/LabNotes/SexChrGenes.txt" 
} else {
  excludedFeaturesFile = NA
}
author <- "Heath O'Brien"                                # author of the statistical analysis/report

workDir <- paste("~/BTSync/FetalRNAseq/Counts", projectName, sep='/')      # working directory for the R session

rawDir <- "~/BTSync/FetalRNAseq/Counts/raw"                                      # path to the directory containing raw counts files

targetFile <- "~/BTSync/FetalRNAseq/LabNotes/MvsFmac.txt"

SequencingCentreFile <- "~/BTSync/FetalRNAseq/LabNotes/SampleProgress.txt"

if (BrainBank == 'HDBRexpression'){
  SampleInfoFile <- "~/BTSync/FetalRNAseq/LabNotes/HDBRsample_info.txt"
} else {
  SampleInfoFile <- "~/BTSync/FetalRNAseq/LabNotes/sample_info.txt"
}
featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")# NULL if no feature to remove

if (BrainBank == 'HDBRexpression'){
  batch <- c()
}else if ( pcw ) {
  batch <- c("Centre", "RIN", "PCW")                # blocking factor: NULL (default) or "batch" for example
} else {
  batch <- c("Centre", "RIN")                # blocking factor: NULL (default) or "batch" for example
}
interact <- 0

varInt <- "Sex"                                    # factor of interest
condRef <- "Female"                                      # reference biological condition

pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

cpmCutoff <- 1                                       # counts-per-million cut-off to filter low counts
gene.selection <- "pairwise"                         # selection of the features in MDSPlot
normalizationMethod <- "TMM"                         # normalization method: "TMM" (default), "RLE" (DESeq) or "upperquartile"

colors <- c("dodgerblue","firebrick1",               # vector of colors of each biological condition on the plots
            "MediumVioletRed","SpringGreen")


################################################################################
###                             running script                               ###
################################################################################
LabNotes <- "~/BTSync/FetalRNAseq/LabNotes/"
dir.create(workDir)
setwd(workDir)
library(devtools)
load_all(pkg = "~/BTSync/Code/R/SARTools")
library(dplyr)
library(tidyr)
library(readr)
source(paste0(LabNotes, "R/FormatGGplot.R"))
source(paste0(LabNotes, 'R/AnalyseDE.R'))

# checking parameters
checkParameters.edgeR(projectName=projectName,author=author,targetFile=targetFile,
                      rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                      condRef=condRef,batch=batch,alpha=alpha,pAdjustMethod=pAdjustMethod,
                      cpmCutoff=cpmCutoff,gene.selection=gene.selection,
                      normalizationMethod=normalizationMethod,colors=colors)

# loading target file
target <- GetTarget(PCW_cutoff[1], PCW_cutoff[2], RIN_cutoff=RIN_cutoff, exclude=exclude, BrainBank=BrainBank)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)
if (! is.na(excludedFeaturesFile)) {
  excludedFeatures <- read_csv(excludedFeaturesFile, col_names = FALSE)
  counts <- counts[!rownames(counts) %in% excludedFeatures$X1, ]
}

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# edgeR analysis
out.edgeR <- run.edgeR(counts=counts, target=target, varInt=varInt, condRef=condRef,
                       batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)

# MDS + clustering
exploreCounts(object=out.edgeR$dge, group=target[,varInt], gene.selection=gene.selection, col=colors)

# summary of the analysis (boxplots, dispersions, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.edgeR(out.edgeR, group=target[,varInt], counts=counts, alpha=alpha, col=colors)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.edgeR(target=target, counts=counts, out.edgeR=out.edgeR, summaryResults=summaryResults,
                  majSequences=majSequences, workDir=workDir, projectName=paste0("MvsF_", PCW_cutoff, collapse='_'), author=author,
                  targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                  condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod, colors=colors,
                  gene.selection=gene.selection, normalizationMethod=normalizationMethod)


################################################################################
###                          additional steps                                ###
################################################################################

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
MalevsFemale.up <- read.delim("tables/MalevsFemale.up.txt", check.names=FALSE)
MalevsFemale.up <- dplyr::select(MalevsFemale.up, Id, Female, Male, FC, log2FoldChange, pvalue, padj)
MalevsFemale.up <- bind_cols(GetGeneIDs(MalevsFemale.up$Id), MalevsFemale.up)
write.table(MalevsFemale.up, file="tables/MaleUp.txt", sep="\t", quote=FALSE, row.names=FALSE)
MalevsFemale.down <- read.delim("tables/MalevsFemale.down.txt", check.names=FALSE)
MalevsFemale.down <- dplyr::select(MalevsFemale.down, Id, Female, Male, FC, log2FoldChange, pvalue, padj)
MalevsFemale.down <- bind_cols(GetGeneIDs(MalevsFemale.down$Id), MalevsFemale.down)
write.table(MalevsFemale.down, file="tables/FemaleUp.txt", sep="\t", quote=FALSE, row.names=FALSE)
DEgenes <- nrow(MalevsFemale.down) + nrow(MalevsFemale.up)
MalevsFemale.complete <- read.delim("tables/MalevsFemale.complete.txt", check.names=FALSE)
filter(MalevsFemale.complete, ! is.na(padj)) %>%
  write_tsv("tables/MalevsFemale.background.txt")

#write summary of analysis to file
summary <- data.frame(Method=c('EdgeR'),
                      BrainBank=c(BrainBank), 
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
                      test=c("None"),
                      model=c(ifelse(length(interact)==0, '+', '*')), 
                      CooksCutoff=c("None"),
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



