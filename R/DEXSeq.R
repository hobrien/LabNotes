library(readr)
library(dplyr)
library(DEXSeq)
library(BiocParallel)
BPPARAM = MulticoreParam(workers=8)

flattenedFile <- "/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gff"

countFiles <- list.files("/c8000xd3/rnaseq-heath/Counts/DEXSeq", pattern=".chr.nsort.dex_counts.txt$", full.names=TRUE)
sampleTable <- read_tsv("/c8000xd3/rnaseq-heath/Counts/DEXSeq/target.txt")
sample_info <- read_tsv("~/LabNotes/sample_info.txt")
sampleTable <- left_join(sampleTable, dplyr::select(sample_info, BrainBankID, Sex, PCW, RIN), by = c("label" = "BrainBankID"))
rownames(sampleTable) <- sampleTable$label
#sampleTable <-sampleTable[1:6,]
sampleTable <- as.data.frame(sampleTable)

dxd = DEXSeqDataSetFromHTSeq(countFiles, #[1:6],
                             sampleData=sampleTable,
                             design= ~ sample + exon + Sex:exon,
                             flattenedfile=flattenedFile
                             )
formulaFullModel    =  ~ sample + exon + centre:exon + PCW:exon + RIN:exon + Sex:exon 
formulaReducedModel    =  ~ sample + exon + centre:exon + PCW:exon + RIN:exon
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM ) #, formula = formulaFullModel, )
png(filename="/c8000xd3/rnaseq-heath/Counts/DEXSeq/DispEsts.png", bg="transparent", width=300, height=300, units="px")
plotDispEsts( dxd )
dev.off()
dxd = testForDEU( dxd, BPPARAM=BPPARAM)
#                  reducedModel = formulaReducedModel,
#                  fullModel = formulaFullModel,
#                  
#)

dxd = estimateExonFoldChanges( dxd, fitExpToVar="Sex", BPPARAM=BPPARAM)
dxr1 = DEXSeqResults( dxd )
write.table(dxr1, file = "/c8000xd3/rnaseq-heath/Counts/DEXSeq/results.txt")
table ( dxr1$padj < 0.1 ) # number of significant features (FDR < 0.1)
table ( tapply( dxr1$padj < 0.05, dxr1$groupID, any ) ) # number of genes affected
png(filename="/c8000xd3/rnaseq-heath/Counts/DEXSeq/plotMA.png", bg="transparent", width=300, height=300, units="px")
plotMA( dxr1, cex=0.8 )
dev.off()
DEXSeqHTML( dxr1, FDR=0.05, color=c("#FF000080", "#0000FF80"), fitExpToVar="Sex", path="/c8000xd3/rnaseq-heath/Counts/DEXSeq", BPPARAM=BPPARAM )
