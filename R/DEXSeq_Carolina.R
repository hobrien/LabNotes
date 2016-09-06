library(readr)
library(dplyr)
library(DEXSeq)
library(BiocParallel)
BPPARAM = MulticoreParam(workers=4)

flattenedFile <- "~/BTSync/FetalRNAseq/Counts/Carolina/CACNA1C.gff"

countFiles <- list.files("~/BTSync/FetalRNAseq/Counts/Carolina/", pattern="chr.dex_counts.txt$", full.names=TRUE)
sampleTable <- read_tsv("~/BTSync/FetalRNAseq/Counts/Carolina/target.txt")
rownames(sampleTable) <- sampleTable$label
#sampleTable <-sampleTable[1:6,]
sampleTable <- as.data.frame(sampleTable)

dxd = DEXSeqDataSetFromHTSeq(countFiles, #[1:6],
                             sampleData=sampleTable,
                             design= ~ sample + exon + region:exon,
                             flattenedfile=flattenedFile
                             )
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM ) #, formula = formulaFullModel, )
png(filename="~/BTSync/FetalRNAseq/Counts/Carolina/DispEsts.png", bg="transparent", width=600, height=600, units="px")
plotDispEsts( dxd )
dev.off()
dxd = testForDEU( dxd, BPPARAM=BPPARAM)
#                  reducedModel = formulaReducedModel,
#                  fullModel = formulaFullModel,
#                  
#)

dxd = estimateExonFoldChanges( dxd, fitExpToVar="region", BPPARAM=BPPARAM)
dxr1 = DEXSeqResults( dxd )
write.table(dxr1, file = "~/BTSync/FetalRNAseq/Counts/Carolina/results.txt")
table ( dxr1$padj < 0.1 ) # number of significant features (FDR < 0.1)
table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) ) # number of genes affected
png(filename="~/BTSync/FetalRNAseq/Counts/Carolina/plotMA.png", bg="transparent", width=600, height=600, units="px")
plotMA( dxr1, cex=0.8 )
dev.off()
DEXSeqHTML( dxr1, FDR=0.1, color=c("#FF000080", "#0000FF80"), fitExpToVar="region", path="~/BTSync/FetalRNAseq/Counts/Carolina/", BPPARAM=BPPARAM )
save.image(file="~/BTSync/FetalRNAseq/Counts/Carolina/results.RData")