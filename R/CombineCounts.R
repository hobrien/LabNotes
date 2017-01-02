library(readr)
library(dplyr)
nameStem <- args[1]
folder <- '/c8000xd3/rnaseq-heath/Mappings/'
AllFiles <- list.files(folder)
Counts <- data.frame('Feature'=c(), 'Count'=c())
for (fileName in AllFiles[grepl(paste0(nameStem, '-'), AllFiles)]) {
  input <- read_delim(paste0(folder, fileName, '/BAM/', fileName, '.chr.counts.txt'), 
             "\t", escape_double = FALSE, col_names = FALSE, 
             trim_ws = TRUE)
  if (length(Counts) == 0) {
    Counts <- rbind(Counts, input)
  } else {
    Counts <- join(Counts, input, by= 'X1')
  }    
}
Counts$Sum <- rowSums(Counts[,-1])
Counts <-Counts[,c(1, ncol(Counts))]
Counts$Sum <- as.integer(Counts$Sum)
write_delim(Counts, 
            paste0(folder, nameStem, '/BAM/', nameStem, '.chr.counts.txt'),
            delim='\t',
            col_names=FALSE
            )
