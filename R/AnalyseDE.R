library(tidyverse)
library(GO.db)
#library(DOSE)
library(org.Hs.eg.db)
#library(topGO)
#library(GSEABase)
#library(clusterProfiler)
source("~/BTSync/FetalRNAseq/LabNotes//R/FormatGGplot.R")

GetFoldChanges <- function() {
  Fold_changes <- read_delim("~/BTSync/FetalRNAseq/Counts/MvsF_12_HDBR_FDR_0.1_edgeR/tables/MalevsFemale.complete.txt",
                             "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    mutate(sig = ifelse(padj < 0.1, 1, 0)) %>%
    dplyr::select(Id, log2FoldChange, sig) %>%
    mutate('week' = 12) %>%
    bind_rows(read_delim("~/BTSync/FetalRNAseq/Counts/MvsF_13_HDBR_excl_16491_FDR_0.1_edgeR/tables/MalevsFemale.complete.txt",
                         "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      mutate(sig = ifelse(padj < 0.1, 1, 0)) %>%
      dplyr::select(Id, log2FoldChange, sig) %>%
      mutate('week' = 13)
    ) %>%
    bind_rows(read_delim("~/BTSync/FetalRNAseq/Counts/MvsF_14_HDBR_FDR_0.1_edgeR/tables/MalevsFemale.complete.txt",
                         "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      mutate(sig = ifelse(padj < 0.1, 1, 0)) %>%
      dplyr::select(Id, log2FoldChange, sig) %>%
      mutate('week' = 14)
    ) %>%
    bind_rows(read_delim("~/BTSync/FetalRNAseq/Counts/MvsF_15_17_HDBR_excl_15641_FDR_0.1_edgeR/tables/MalevsFemale.complete.txt",
                         "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      mutate(sig = ifelse(padj < 0.1, 1, 0)) %>%
      dplyr::select(Id, log2FoldChange, sig) %>%
      mutate('week' = 15.5)
    ) %>%
      bind_rows(read_delim("~/BTSync/FetalRNAseq/Counts/MvsF_17_20_HDBR_excl_18432_FDR_0.1_edgeR/tables/MalevsFemale.complete.txt",
                           "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      mutate(sig = ifelse(padj < 0.1, 1, 0)) %>%
      dplyr::select(Id, log2FoldChange, sig) %>%
      mutate('week' = 18)
    ) %>%
    mutate(Id = gsub('\\..*', '', Id))
  Fold_changes
}

PlotFC <- function(genes, set_name) {
  plot <- ggplot(genes, aes(x=week, y=log2FoldChange, group=Id)) +
    geom_line(alpha=.1) +
    geom_point(colour='red', size=1, alpha=.25, data=filter(genes, sig==1)) +
    scale_x_continuous(breaks=c(12,13,14,15.5,18), labels=c('12', '13', '14', '15-16', '17-19')) +
    tufte_theme() +
    ylab("<--- Female-biased        Male-biased --->") +
    xlab("Post Conception Week") +
    theme(axis.title.y=element_text(size=9)) +
    ggtitle(paste("Log fold differences in", set_name))
  plot        
}

PlotFC_by_ChrType <- function(genes, set_name) {
  plot <- ggplot(genes, aes(x=week, y=log2FoldChange, group=Id)) +
    geom_line(alpha=.25, aes(colour=ChrType)) +
    geom_point(colour='red', size=1, alpha=.25, data=filter(genes, sig==1)) +
    scale_x_continuous(breaks=c(12,13,14,15.5,18), labels=c('12', '13', '14', '15-16', '17-19')) +
    scale_colour_manual(values=c('black', 'orange', 'blue')) +
    tufte_theme() +
    ylab("<--- Female-biased        Male-biased --->") +
    xlab("Post Conception Week") +
    theme(axis.title.y=element_text(size=9)) +
    ggtitle(paste("Log fold differences in", set_name))
  plot        
}
min=12
max=20

GetTarget <- function(min, max, RIN_cutoff=0, exclude=NA, BrainBank='HDBR', varInt='Sex') {
  PCW_cutoff <- c(min, max)
  targetFile <- "~/BTSync/FetalRNAseq/LabNotes/MvsFmac.txt"
  
  SequencingCentreFile <- "~/BTSync/FetalRNAseq/LabNotes/SampleProgress.txt"
  
  if (BrainBank == 'HDBRexpression'){
    SampleInfoFile <- "~/BTSync/FetalRNAseq/LabNotes/HDBRsample_info.txt"
  } else {
    SampleInfoFile <- "~/BTSync/FetalRNAseq/LabNotes/sample_info.txt"
  }
  
  target <- read_tsv(targetFile)                        # path to the design/target file
  
  sample_info <- read_tsv(SampleInfoFile, col_types=cols(BrainBankID='c', Sex='c', PCW='n', RIN='n'))
  if ( BrainBank == 'HDBRexpression') {
    target <- right_join(target, sample_info, by = c("label" = "BrainBankID"))
  } else {
    target <- left_join(target, dplyr::select(sample_info, BrainBankID, Sex, PCW, RIN), by = c("label" = "BrainBankID"))
    SequencingCentre <- read.delim(SequencingCentreFile)
    target <- left_join(target, dplyr::select(SequencingCentre, sample, Centre), by = c("label" = "sample"))
  }
  target <- arrange(target, Sex)
  if (!is.na(RIN_cutoff)) {
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
  
  # covert PCW+0 samples to PCW-1
  if (BrainBank != 'HDBRexpression') { 
    newer_samples = c('11875', 
                      '12107', 
                      '12545',
                      '12546',
                      '12994',
                      '13142',
                      '12993',
                      '13008') #these are already the correct week
    #target$PCW <- ifelse(target$PCW %% 1 == 0 & ! target$label %in% newer_samples,
    #                     target$PCW -1, floor(target$PCW))
  }
  target$PCW <- floor(target$PCW)
  target <- droplevels(target)
  target <- mutate(target, Sex = factor(ifelse(Sex == 'unknown', NA, Sex)))
  target <- as.data.frame(target)
  target[,varInt] <- as.factor(target[,varInt])
  target
}

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

GetCounts <- function(){
  MalevsFemale_complete <- read_delim("~/BTSync/FetalRNAseq/Counts/MvsF_12_20_HDBR_excl_15641_18432_16491_FDR_0.1_edgeR/tables/MalevsFemale.complete.txt",
                                      "\t", escape_double = FALSE, trim_ws = TRUE)
  SexChrGenes <- read_delim("/Users/heo3/BTSync/FetalRNAseq/LabNotes/SexChrGenes.txt", 
                            delim='\t', 
                            col_names=c('Id', 'Chr'), 
                            col_types=cols(Id='c', Chr='c')
  )

  MalevsFemale_complete <- MalevsFemale_complete %>%
    full_join(SexChrGenes) %>%
    mutate(ChrType = ifelse(is.na(Chr), 'autosomal', Chr), Id = gsub('\\..*', '', Id))
  MalevsFemale_complete
}

GetGeneSets <- function() {
  gene_sets <- read_delim("~/BTSync/FetalRNAseq/GSEA/MSigDB/PollenEtAl.txt", 
                          "\t", escape_double = FALSE, col_names=c('gene_name', 'set', 'reference'), 
                          trim_ws = TRUE) %>%
    bind_rows(read_delim("~/BTSync/FetalRNAseq/GSEA/MSigDB/DarmanisEtAl.txt", 
                         "\t", escape_double = FALSE, col_names=c('gene_name', 'set', 'reference'), 
                         trim_ws = TRUE))
  gene_sets <- full_join(gene_sets, 
                         bitr(gene_sets$gene_name, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db"),
                         by=c('gene_name' = 'SYMBOL')
  )
  BG <- MalevsFemale_complete %>%
    filter(! is.na(FC)) %>%
    dplyr::select(Id) %>%
    mutate(set="Background", reference = NA)

  gene_sets <- bitr(BG$Id, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db") %>%
    full_join(BG, by=c('ENSEMBL' = 'Id')) %>%
    mutate(gene_name = SYMBOL) %>%
    dplyr::select(gene_name, set, reference, ENSEMBL) %>%
    bind_rows(gene_sets)
  gene_sets
}

PlotExpression<-function(geneID, fileName, id="A value that hopefully isn't in the dataset") {
  Ensembl_Id <- bitr(geneID, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")[,2]
  data <- filter(fileName, Id == Ensembl_Id) %>%  
    dplyr::select(starts_with('norm')) %>%
    gather() %>%
    separate(key, into=c('norm', 'label'), sep='[.]') %>%
    dplyr::select(label, value) %>%
    left_join(target)
  plot<-  ggplot(subset(data, label != id), aes(x=PCW, y=value, colour=Sex)) + 
    geom_jitter(height = 0, width=.1, alpha=.75) + 
    geom_point(data=subset(data, label==id), colour='orange') +
    geom_smooth() +
    ylab("normalised counts") +
    tufte_theme() +
    scale_colour_brewer(type = "qual", palette = 6) +
    ggtitle(geneID) +
    theme(legend.position=c(0.1,.9)) +
    theme(plot.background=element_blank())
  plot
}

PlotTimepoint<-function(geneID, fileName) {
  data <- fileName %>%  
    dplyr::select(starts_with('norm')) %>%
    gather() %>%
    separate(key, into=c('norm', 'label'), sep='[.]') %>%
    dplyr::select(label, value) %>%
    left_join(target)
  mean <- fileName %>% 
    dplyr::select(Male, Female) %>%
    gather()
  plot<-  ggplot(data, aes(x=Sex, y=value, colour=Sex)) + 
    geom_errorbar(aes(x=key, ymin=value, ymax=value), colour='black', size=1, width=.5, data=mean) +
    geom_jitter(height = 0, width=.1, alpha=.75) + 
    ylab("normalised counts") +
    xlab('') +
    tufte_theme() +
    scale_colour_brewer(type = "qual", palette = 6) +
    ggtitle(geneID) +
    theme(plot.background=element_blank())
  plot
}
