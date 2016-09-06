library(JunctionSeq)

decoder <- read.table("~/Downloads/JctSeqData/inst/extdata/annoFiles/decoder.bySample.txt",
                      header=TRUE,
                      stringsAsFactors=FALSE);
gff.file <- "~/Downloads/JctSeqData/inst/extdata/cts/withNovel.forJunctionSeq.gff.gz"
countFiles.noNovel<-"~/Downloads/JctSeqData/inst/extdata/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz"
countFiles<- "~/Downloads/JctSeqData/inst/extdata/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"
gff.file<-"~/Downloads/JctSeqData/inst/extdata/tiny/withNovel.forJunctionSeq.gff.gz"
countFiles<- paste0("~/Downloads/JctSeqData/inst/extdata/tiny/",
                    decoder$sample.ID,
                    "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

decoder <- read.table("~/BTSync/FetalRNAseq/Counts/Carolina/decoder.bySample.txt",
                      header=TRUE,
                      stringsAsFactors=FALSE);
gff.file <- "~/BTSync/FetalRNAseq/Counts/Carolina/withNovel.forJunctionSeq.gff.gz"
countFiles<- paste0("~/BTSync/FetalRNAseq/Counts/Carolina/",
                    decoder$sample.ID,
                    "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                               sample.names = decoder$sample.ID,
                               condition=factor(decoder$group.ID),
                               flat.gff.file = gff.file,
                               nCores = 4,
                               analysis.type = "junctionsAndExons"
)

buildAllPlots(jscs=jscs,
              gene.list = c("ENSG00000151067.20", 
                            "ENSG00000256837.1", 
                            "ENSG00000256257.1", 
                            "ENSG00000256025.1", 
                            "ENSG00000256721.1", 
                            "ENSG00000256769.1",
                            "ENSG00000246627.6"),
              outfile.prefix = "~/BTSync/FetalRNAseq/Counts/Carolina/plotsJCT-2/",
              use.plotting.device = "png",
              FDR.threshold = 0.01)

buildAllPlotsForGene(jscs=jscs,
                     geneID = "ENSG00000151067.20",
              outfile.prefix = "~/BTSync/FetalRNAseq/Counts/Carolina/plotsJCT_genes/",
              use.plotting.device = "png",
              expr.plot = FALSE, normCounts.plot = FALSE,
              rExpr.plot = TRUE, rawCounts.plot = FALSE,
              without.TX=FALSE,
              colorList = list(SIG.FEATURE.FILL.COLOR = "red"),
              colorRed.FDR.threshold = 0.001)
