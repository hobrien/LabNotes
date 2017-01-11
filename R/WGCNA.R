################################################################################
### R script to run WCGNA analyses, Inspired by SARTools (Hugo Varet)
### Heath O'Brien
### 10 Jan 2017
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
rm(list=ls())                                        # remove all the objects from the R session

workDir <- "~/BTSync/FetalRNAseq/WGCNA"      # working directory for the R session

projectName <- "MvsF_12_20_Cooks.75_excl_15641_18432_16491_FDR.1_new_VST1_comBat_QN"                         # name of the project
author <- "Heath O'Brien"                                # author of the statistical analysis/report

workDir <- paste("~/BTSync/FetalRNAseq/WGCNA", projectName, sep='/')      # working directory for the R session

inputDir <- "~/BTSync/FetalRNAseq/Counts/MvsF_12_20_Cooks.75_excl_15641_18432_16491_FDR.1_new/tables"

#inputCounts <- "VST_filtered.txt" #unfiltered Variance Stabalising Transformation of normalised counts
inputCounts <- "VST_filtered.txt" #filter out genes with counts < 10 in more than 10% of samples
#inputCounts <- "VST_filtered2.txt" #filter out genes with counts < 10 in more than 50% of samples

batch <- c("Centre", "RIN")                  
quantile_norm <- TRUE

################################################################################
###                             running script                               ###
################################################################################

library(WGCNA)
library(readr)
library(dplyr)
library(sva)
library(preprocessCore)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

dir.create(workDir)
dir.create(paste(workDir, 'figures', sep='/'))
setwd(workDir)
getwd()

# Prepare data for analysis
MalevsFemale <- read_delim(paste(inputDir, inputCounts, sep='/'), "\t", escape_double = FALSE, trim_ws = TRUE)
#move Id to rownames
rownames(MalevsFemale) <- MalevsFemale$Id
MalevsFemale <- select(MalevsFemale, -1) 

dim(MalevsFemale)
names(MalevsFemale)

#Add trait info
sample_info <- read.delim("~/BTSync/FetalRNAseq/LabNotes/sample_info.txt") %>%
  select(BrainBankID, Sex, PCW, RIN)
sample_info <- read.delim("~/BTSync/FetalRNAseq/LabNotes/SampleProgress.txt") %>%
  select(sample, Centre) %>%
  left_join(sample_info, by = c("sample" = "BrainBankID"))
SampleNames <- colnames(MalevsFemale);
traitRows = match(SampleNames, sample_info$sample);
datTraits = sample_info[traitRows, -1];
rownames(datTraits) = sample_info[traitRows, 1];
collectGarbage()

datTraits$Sex <- as.numeric(as.factor(datTraits$Sex))
datTraits$Centre <- as.numeric(as.factor(datTraits$Centre))

#Correct for batch effects (see http://jtleek.com/genstats/inst/doc/02_13_batch-effects.html)
if ('Centre' %in% batch) {
    batch <- datTraits$Centre
    modcombat = model.matrix(~1, data=datTraits)
    MalevsFemale = ComBat(dat=MalevsFemale, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
}

if (quantile_norm) {
    MalevsFemale <- normalize.quantiles(as.matrix(MalevsFemale),copy=TRUE)
}

#transpose, remove metadata about genes
datExpr0 <- as.data.frame(t(MalevsFemale))
rownames(datExpr0) <- colnames(MalevsFemale)
colnames(datExpr0) <- rownames(MalevsFemale)


#test for gene and sample missingness
gsg = goodSamplesGenes(datExpr0, verbose = 3);
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)



traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
png(filename="figures/clusterSamples.png",width=1800,height=1800,res=300) 
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:

# Scale-free topology fit index as a function of the soft-thresholding power
png(filename="figures/softThreshold.png",width=1800,height=1800,res=300) 
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

net = blockwiseModules(datExpr0, power = 10, maxBlockSize = 20000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "MalevsFemaleFilter_8TOM",
                       verbose = 3)

table(net$colors)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
png(filename="figures/modules.png",width=1800,height=1800,res=300) 
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

png(filename="figures/heatmap.png",width=1800,height=1800,res=300) 
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

save.image(file=paste0(projectName, ".RData"))

