library(tidyverse)

# Oh yes, thought of a catchy (and hopefully memorable) name for the dataset: 
# GENEX-FB (for GENe EXpression in the Fetal Brain). 
# This dataset would be GENEX-FB1: Sex biases. 
# The larger dataset for eQTL / TWAS analysis (I reckon we'll get up to 150) would be GENEX-FB2: 
# Genotypic effects. (If we grew the sample, we'd call it FB3 etc). 
# We could also use GENEX for the adult brain samples (E.g. GENEX-AC (adult caudate)!

read_delim("~/BTSync/FetalRNAseq/LabNotes/Results/BG12_19.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(Id, SYMBOL, starts_with('norm')) %>% write_tsv("~/BTSync/FetalRNAseq/LabNotes/R/GENEX-FB1/Data/counts.txt")

read_delim("~/BTSync/FetalRNAseq/LabNotes/Results/BG12_19.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(Id, SYMBOL, Male, Female, log2FoldChange, pvalue, padj) %>% 
  mutate(ageBin='12-19') %>%
  bind_rows(
    read_delim("~/BTSync/FetalRNAseq/LabNotes/Results/BG12.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, SYMBOL, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='12')
  ) %>%
  bind_rows(
    read_delim("~/BTSync/FetalRNAseq/LabNotes/Results/BG13.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, SYMBOL, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='13')
  ) %>%
  bind_rows(
    read_delim("~/BTSync/FetalRNAseq/LabNotes/Results/BG14.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, SYMBOL, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='14')
  ) %>%
  bind_rows(
    read_delim("~/BTSync/FetalRNAseq/LabNotes/Results/BG15_16.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, SYMBOL, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='15-16')
  ) %>%
  bind_rows(
    read_delim("~/BTSync/FetalRNAseq/LabNotes/Results/BG17_19.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, SYMBOL, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='17-19')
  ) %>%
  write_tsv("~/BTSync/FetalRNAseq/LabNotes/R/GENEX-FB1/Data/fitted.txt")
