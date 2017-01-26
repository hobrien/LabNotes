library(readr)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
Counts <- data.frame('Feature'=c(), 'Count'=c())
for (fileName in args) {
  input <- read_delim(paste0(folder, fileName, '/BAM/', fileName, '.chr.counts.txt'), 
             "\t", escape_double = FALSE, col_names = FALSE, 
             trim_ws = TRUE)
  if (length(Counts) == 0) {
    Counts <- rbind(Counts, input)
  } else {
    Counts <- full_join(Counts, input, by= 'X1')
  }    
}
Counts$Sum <- rowSums(Counts[,-1])
Counts <-Counts[,c(1, ncol(Counts))]
Counts$Sum <- as.integer(Counts$Sum)
write_delim(Counts, 
            sub('.*/', '~/Counts/', args[1]),
            delim='\t',
            col_names=FALSE
            )
