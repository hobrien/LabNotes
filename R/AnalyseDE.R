library(tidyverse)
library(GO.db)
#library(DOSE)
library(org.Hs.eg.db)
#library(topGO)
#library(GSEABase)
#library(clusterProfiler)
source("~/BTSync/FetalRNAseq/LabNotes//R/FormatGGplot.R")

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

PlotFC <- function(genes, set_name) {
  plot <- ggplot(genes, aes(x=week, y=log2FoldChange, group=gene_name)) +
    geom_line(alpha=.1) +
    geom_point(colour='red', size=1, alpha=.25, data=filter(genes, sig==1)) +
    scale_x_continuous(breaks=c(12,13,14,15.5,18), labels=c('12', '13', '14', '15-16', '17-19')) +
    tufte_theme() +
    ylab("<--- Female-biased        Male-biased --->") +
    theme(axis.title.y=element_text(size=9)) +
    ggtitle(paste("Log fold differences in", set_name, "genes"))
  plot        
}

PlotFC_by_ChrType <- function(genes, set_name) {
  plot <- ggplot(genes, aes(x=week, y=log2FoldChange, group=gene_name)) +
    geom_line(alpha=.25, aes(colour=ChrType)) +
    geom_point(colour='red', size=1, alpha=.25, data=filter(genes, sig==1)) +
    scale_x_continuous(breaks=c(12,13,14,15.5,18), labels=c('12', '13', '14', '15-16', '17-19')) +
    scale_colour_manual(values=c('black', 'orange', 'blue')) +
    tufte_theme() +
    ylab("<--- Female-biased        Male-biased --->") +
    theme(axis.title.y=element_text(size=9)) +
    ggtitle(paste("Log fold differences in", set_name, "genes"))
  plot        
}
