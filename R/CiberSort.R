library(tiyverse)
library(clusterProfiler)
library(ggdendro)
library(gridExtra)
LabNotes <- "~/BTSync/FetalRNAseq/LabNotes/"
source(paste0(LabNotes, 'R/AnalyseDE.R'))
Counts <- GetCounts() %>% dplyr::select(Id, starts_with('norm.'))
Counts <- subset(Counts, ! is.na(Counts[,2]))
matrix_pure <- read_delim(paste0(LabNotes, "Cibersort/matrix_pure.txt"),
"\t", escape_double = FALSE, trim_ws = TRUE)
Counts2 <- bitr(matrix_pure$`!Sample_title`, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db") %>%
  inner_join(Counts, by=c('ENSEMBL' = 'Id')) %>% dplyr::select(-ENSEMBL)
colnames(Counts2)[1] <- "!Sample_title"
write_tsv(Counts2, paste0(LabNotes, "Cibersort/matrix_mix.txt"))


full_matrix <- read_delim(paste0(LabNotes, "Cibersort/Pollen500.txt"),
                          "\t", escape_double = FALSE, trim_ws = TRUE) %>% as.data.frame()
rownames(full_matrix) <- full_matrix$`!Sample_title`
full_matrix <- dplyr::select(full_matrix, -`!Sample_title`)
heatmap(t(as.matrix(full_matrix)))
class <- data.frame(sampleID = colnames(full_matrix),
                    class=c(rep(NPC, 15), )
cluster <- hclust(dist(t(as.matrix(full_matrix))))
ggdendrogram(cluster)+geom_hline(yintercept =131)
groups <-as.data.frame(cutree(cluster, h=131)) 
colnames(groups) <-'cluster'

groups$NPC<-ifelse(groups$cluster == 1, 1, 2)
groups$Radial_glia<-ifelse(groups$cluster == 2, 1, 2)
groups$New_neuron<-ifelse(groups$cluster == 3, 1, 2)
groups$Mature_neuron<-ifelse(groups$cluster == 4, 1, 2)
t(dplyr::select(groups, -cluster)) %>% as.data.frame() %>% write.table(paste0(LabNotes, "Cibersort/class_full.txt"), sep='\t', col.names=F, row.names=T, quote=F)

cell_types <- read_delim(paste0(LabNotes, "Cibersort/CIBERSORT.Output_Job6.txt"),
                          "\t", escape_double = FALSE, trim_ws = TRUE) %>% separate(`Input Sample`, c("x", "Sample")) %>% dplyr::select(-x)
cell_types <- cell_types %>% gather(CellType, Proportion, -Sample, -`P-value`, -`Pearson Correlation`, -`RMSE`)

target <-GetTarget(12, 20)
cell_types <- inner_join(cell_types, target, by=c("Sample" = "label"))

PlotCellTypes <- function(cell_types, MorF, Min, Max) {
  p1 <- filter(cell_types, Sex == MorF, PCW >= Min, PCW < Max) %>% 
  ggplot(aes(x=Sample, y=Proportion, fill=CellType))+geom_bar(stat="identity") + theme(legend.position="none")
  p1
}
layout(matrix(seq(1,10), 2))
grid.arrange(PlotCellTypes(cell_types, 'Male', 12, 13),
             PlotCellTypes(cell_types, 'Male', 13, 14),
             PlotCellTypes(cell_types, 'Male', 14, 15),
             PlotCellTypes(cell_types, 'Male', 15, 17),
             PlotCellTypes(cell_types, 'Male', 17, 20),
             PlotCellTypes(cell_types, 'Female', 12, 13),
             PlotCellTypes(cell_types, 'Female', 13, 14),
             PlotCellTypes(cell_types, 'Female', 14, 15),
             PlotCellTypes(cell_types, 'Female', 15, 17),
             PlotCellTypes(cell_types, 'Female', 17, 20),
             ncol=5)


