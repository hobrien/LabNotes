library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(tufte)
source("~/BTSync/FetalRNAseq/LabNotes/R/FormatGGplot.R")

GetCounts<-function(count_file, riskVar_file, ExpressedSNP_file) {
  counts <- read_tsv(count_file, 
                     col_names=c('chr', 'pos', 'ref_base', 'ref', 'alt', 'sample'),
                     col_types = cols(chr="c", 
                                      pos="i", 
                                      ref_base="c", 
                                      ref="i",
                                      alt="i",
                                      sample="c"
                     )
  )
  riskVar <- read_tsv(riskVar_file, 
                      col_names=c('chr', 'pos', 'riskVar', 'sample'), 
                      col_types =cols(chr='c', 
                                      pos='i', 
                                      riskVar='c', 
                                      sample='c'
                      )
  )
  ExpressedSNP <- read_tsv(ExpressedSNP_file,
                           col_names=c('chr', 'pos', 'ExpressedSNP', 'sample'),
                           col_types =cols(chr='c', 
                                           pos='i', 
                                           ExpressedSNP='c', 
                                           sample='c'
                           )
  )
  counts <- left_join(counts, riskVar[,c(3,4)]) %>% 
    left_join(ExpressedSNP[,c(3,4)])
}

PlotDistortion <- function(counts) {
  p1 <- ggplot(counts, aes(x=riskVar, y=alt/(ref+alt), size=ref+alt)) +
    geom_jitter(alpha=.2, height=0, width=.25) +
    geom_abline(intercept=0.5, slope=0, colour='red') +
    scale_y_continuous(limits=c(-.05, 1.05), breaks=c(0,.2,.4,.6,.8,1)) +
    ylab('% Alternate Allele') +
    xlab('GWAS Risk Variant') +
    tufte_theme() +
    theme(legend.position="top")
  print(p1)
}

PlotRatio <- function(counts) {
  p1 <- ggplot(counts, aes(x=riskVar, y=log2(alt/ref), size=ref+alt)) +
    geom_jitter(alpha=.2, height=0, width=.25) +
    geom_abline(intercept=0, slope=0, colour='red') +
    scale_y_continuous(limits=c(-2,2), breaks=c(-2,-1,0,1,2), labels=c(2^-2,2^-1,2^0,2^1,2^2)) +
    ylab('alt/ref ratio') +
    xlab('GWAS Risk Variant') +
    tufte_theme() +
    theme(legend.position="top")
  print(p1)
}

PlotASE <- function(counts) {
  counts <- filter(counts, !is.na(ExpressedSNP))
  p1 <- ggplot(counts, aes(x=riskVar, y=alt/(ref+alt), size=ref+alt)) +
    geom_jitter(alpha=.2, height=0, width=.25) +
    facet_grid(ExpressedSNP ~ .) +
    scale_y_continuous(limits=c(-.1, 1.1), breaks=c(0,.2,.4,.6,.8,1)) +
    ylab('% Alternate Allele') +
    xlab('GWAS Risk Variant') +
    tufte_theme() +
    theme(legend.position="top")
  print(p1)
}

GetAllCounts <- function(snp_pos, riskVar_file, ExpressedSNP_file) {
  count_file <- paste0("~/BTSync/FetalRNAseq/Counts/AlleleSpecific/", snp_pos, ".counts.txt")
  RawCounts<-GetCounts(count_file, riskVar_file, ExpressedSNP_file)
  count_file <- paste0("~/BTSync/FetalRNAseq/Counts/AlleleSpecific/", snp_pos, ".dedup.counts.txt")     
  Deduplicated<-GetCounts(count_file, riskVar_file, ExpressedSNP_file)
  count_file <- paste0("~/BTSync/FetalRNAseq/Counts/AlleleSpecific/", snp_pos, ".clip.counts.txt")
  Clipped<-GetCounts(count_file, riskVar_file, ExpressedSNP_file)
  CombineCounts(RawCounts, Deduplicated, Clipped)
}

CombineCounts <- function(RawCounts, Deduplicated, Clipped) {
   rename(RawCounts, ref.raw = ref, alt.raw = alt) %>%
          full_join(select(Deduplicated, 
                           sample, 
                           ref.dedup = ref, 
                           alt.dedup = alt
                           ), by=c("sample")) %>%
          full_join(select(Clipped, 
                          sample, 
                          ref.clip = ref, 
                          alt.clip = alt
                          ), by=c("sample")) 
}  

PlotCounts <- function(Combined) {
  p1 <- mutate(Combined, raw = ref.raw + alt.raw, 
         dedup = ref.dedup + alt.dedup, 
         clip =  ref.clip + alt.clip) %>%
    select(raw, dedup, clip) %>%
    ggpairs() +tufte_theme()
  print(p1)
}
