#library(tidyverse)
library(dplyr)
library(readr)
library(tidyr)
#LabNotes="/Users/heo3/BTSync/FetalRNAseq/LabNotes/"
LabNotes="~/LabNotes/"
SeqInfo <- read_delim(paste0(LabNotes, "sequences.txt"), 
                     delim='\t',
                     col_names=c('read_file', 'read_group', 'centre', 'folder'),
                     col_types=cols(read_file='c', read_group='c', centre='c', folder='c')
                     )

 
counts_files <- read_delim(paste0(LabNotes, "chr_bam_R1.txt"),
"\t", escape_double = FALSE, col_names = c('file_path', 'counts_file', 'read_path', 'read_file'),
trim_ws = TRUE)

in_stats<- read_delim(paste0(LabNotes, "in_stats_files.txt"), "\t", escape_double = FALSE, col_names = c('file_path', 'in_stats'),
trim_ws = TRUE)

ex_stats<- read_delim(paste0(LabNotes, "stats_files.txt"), "\t", escape_double = FALSE, col_names = c('file_path', 'ex_stats'),
trim_ws = TRUE)

SeqInfo <- filter(SeqInfo, grepl('_1$', read_file) | grepl('_R1_', read_file)) %>%
  inner_join(counts_files) %>% 
  inner_join(in_stats) %>%
  inner_join(ex_stats) 

GetStats <- function(path, in_stats, ex_stats) {
  temp <- read.delim(paste0(path, ex_stats), 
                     header=FALSE, stringsAsFactors=FALSE, skip=5, sep=':'
  )
  temp[4,2] <-strsplit(temp[4,1], ' +')[[1]][4]
  temp[4,1] <- 'Non primary hits'
  temp <- temp[c(6,7,14),]
  temp2 <- read.delim(paste0(path, in_stats), 
                      header=FALSE, stringsAsFactors=FALSE, skip=5, sep=':'
  )
  temp2[4,2] <-strsplit(temp[4,1], ' +')[[1]][4]
  temp2[4,1] <- 'Non primary hits'
  temp <-rbind(temp, c("rDNA",sum(as.numeric(temp2[c(6,7),2]))))
  temp[,1] <- c("Multimapped", "Unique", "Paired", "rDNA")
  temp$V2 <- as.integer(as.character(temp$V2))
  spread(temp, V1, V2)
}
SeqInfo %>%
  write_delim("~/Results/SeqInfo.txt")
Stats <- SeqInfo %>%
  filter(! is.na(file_path)) %>%
  group_by(file_path) %>%
  slice(1L) %>% ungroup() %>% group_by(file_path, in_stats, ex_stats) %>%
  do(GetStats(.$file_path, .$in_stats, .$ex_stats))
write_delim(Stats, "~/Results/Mapping_stats.txt")
