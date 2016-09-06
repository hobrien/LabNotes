library(JunctionSeq)

decoder <- read.table("~/LabNotes/Carolina/decoder.bySample.txt",
                      header=TRUE,
                      stringsAsFactors=FALSE);
gff.file <- "/c8000xd3/rnaseq-heath/Carolina/Counts/withNovel.forJunctionSeq.gff.gz"
countFiles<- paste0("/c8000xd3/rnaseq-heath/Carolina/Counts/",
                    decoder$sample.ID,
                    "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                               sample.names = decoder$sample.ID,
                               condition=factor(decoder$group.ID),
                               flat.gff.file = gff.file,
                               nCores = 8,
                               analysis.type = "junctionsAndExons"
)

buildAllPlots(jscs=jscs,
              outfile.prefix = "/c8000xd3/rnaseq-heath/Carolina/plotsJCT-2/",
              use.plotting.device = "png",
              FDR.threshold = 0.01)

buildAllPlots(jscs=jscs,
              outfile.prefix = "/c8000xd3/rnaseq-heath/Carolina/plotsJCT-3/",
              use.plotting.device = "png",
              FDR.threshold = 0.001)
save.image(file="/c8000xd3/rnaseq-heath/Carolina/junctionSeq.RData")