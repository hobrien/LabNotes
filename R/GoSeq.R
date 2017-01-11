library(goseq)
library(tidyverse)
library(GO.db)
Background <- read_delim("~/BTSync/FetalRNAseq/Counts/MvsF_14_20_noA_Cooks.75_excl_15641_18432_FDR.1_new/tables/Background.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
MaleUp <-     read_delim("~/BTSync/FetalRNAseq/Counts/MvsF_14_20_noA_Cooks.75_excl_15641_18432_FDR.1_new/tables/MaleUp.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
MaleUp <- separate(MaleUp, Id, c('Id', 'Version'), sep="\\.")
assayed.genes <- Background$Id
de.genes <- MaleUp$Id
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes
head(gene.vector)
pwf=nullp(gene.vector,"hg19","ensGene")
GO.wall=goseq(pwf,"hg19","ensGene")
head(GO.wall)
GO.samp=goseq(pwf,"hg19","ensGene",method="Sampling",repcnt=1000)
plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.samp[,1],GO.wall[,1]),2]),
           xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
           xlim=c(-3,0))
abline(0,1,col=3,lty=2)

GO.nobias=goseq(pwf,"hg19","ensGene",method="Hypergeometric")
plot(log10(GO.wall[,2]), log10(GO.nobias[match(GO.nobias[,1],GO.wall[,1]),2]),
           xlab="log10(Wallenius p-values)", ylab="log10(Hypergeometric p-values)",
           xlim=c(-10,0), ylim=c(-10,0))
abline(0,1,col=3,lty=2)

enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,
                                                   method="BH")<.05]
head(enriched.GO)
for(go in enriched.GO[1:10]){
           print(GOTERM[[go]])
           cat("--------------------------------------\n")
}

# I can't be assed to calculate the total counts, so I'm just going to give the sum of fitted males counts plus fitted femal counts. It shouldn't matter
countbias <- Background$Male + Background$Female
pwf.counts=nullp(gene.vector,bias.data=countbias)
GO.counts=goseq(pwf.counts,"hg19","ensGene")
head(GO.counts)
plot(log10(GO.wall[,2]), log10(GO.counts[match(GO.counts[,1],GO.wall[,1]),2]),
     xlab="log10(Length-corrected p-values)", ylab="log10(count-corrected)",
     xlim=c(-12,0), ylim=c(-12,0))
abline(0,1,col=3,lty=2)
enriched.GO=GO.counts$category[p.adjust(GO.counts$over_represented_pvalue,
                                      method="BH")<.05]
head(enriched.GO)
for(go in enriched.GO[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}



FemaleUp <-     read_delim("~/BTSync/FetalRNAseq/Counts/MvsF_14_20_noA_Cooks.75_excl_15641_18432_FDR.1_new/tables/FemaleUp.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
FemaleUp <- separate(FemaleUp, Id, c('Id', 'Version'), sep="\\.")
assayed.genes <- Background$Id
de.genes <- FemaleUp$Id
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes
head(gene.vector)
pwf=nullp(gene.vector,"hg19","ensGene")
GO.wall=goseq(pwf,"hg19","ensGene")
head(GO.wall)

GO.nobias=goseq(pwf,"hg19","ensGene",method="Hypergeometric")
plot(log10(GO.wall[,2]), log10(GO.nobias[match(GO.nobias[,1],GO.wall[,1]),2]),
     xlab="log10(Wallenius p-values)", ylab="log10(Hypergeometric p-values)",
     xlim=c(-10,0), ylim=c(-10,0))
abline(0,1,col=3,lty=2)

enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,
                                      method="BH")<.05]
head(enriched.GO)
for(go in enriched.GO[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}


