library(WGCNA)
library(readr)
library(dplyr)
library(sva)
library(preprocessCore)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

rm(list=ls())

setwd("~/BTSync/FetalRNAseq/WGCNA")
getwd()

# Prepare data for analysis
MalevsFemale <- read_delim("~/BTSync/FetalRNAseq/Counts/MvsF_12_20_noA/tables/VST_filtered.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#move Id to rownames
rownames(MalevsFemale) <- MalevsFemale$Id
MalevsFemale <- select(MalevsFemale, -1) 

dim(MalevsFemale)
names(MalevsFemale)

#Correct for batch effects
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

batch <- datTraits$Centre
mod = model.matrix(~as.factor(Sex) + RIN + PCW + as.factor(batch),data=datTraits)
mod2 = model.matrix(~as.factor(Sex) + RIN + PCW,data=datTraits)
fit = lm.fit(mod,t(MalevsFemale))
hist(fit$coefficients[2,],col=2,breaks=100)
table(datTraits$Centre, datTraits$Sex)
fit2 = lm.fit(mod2,t(MalevsFemale))
hist(fit2$coefficients[2,],col=2,breaks=100)
modcombat = model.matrix(~1, data=datTraits)
MalevsFemale_combat = ComBat(dat=MalevsFemale, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
combat_fit = lm.fit(mod2,t(MalevsFemale_combat))
hist(combat_fit$coefficients[2,],col=2,breaks=100)
plot(fit2$coefficients[2,],combat_fit$coefficients[2,],col=2,
     xlab="Linear Model",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)
qqplot(fit2$coefficients)
MalevsFemale_norm <- normalize.quantiles(as.matrix(MalevsFemale_combat),copy=TRUE)

#transpose, remove metadata about genes
datExpr0 <- as.data.frame(t(MalevsFemale_norm))
rownames(datExpr0) <- colnames(MalevsFemale_combat)
colnames(datExpr0) <- rownames(MalevsFemale_combat)

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


#Add trait info
datTraits$Sex <- as.numeric(as.factor(datTraits$Sex))
datTraits$Centre <- as.numeric(as.factor(datTraits$Centre))

traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

save(datExpr0, file = "MvsF_12_20_noA_filtered-dataInput.RData")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
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

net = blockwiseModules(datExpr0, power = 8, maxBlockSize = 20000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "MalevsFemaleFilter_8TOM",
                       verbose = 3)

table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "MvsF_12_20_noA_filter2_8-auto.RData")

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
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

