library(MatrixEQTL)
library(tidyverse)
library(stringr)
LabNotes="~/BTSync/FetalRNAseq/LabNotes/"
source(paste0(LabNotes, 'R/AnalyseDE.R'))

# Import normalised counts, filter, remove 'norm.', sort
counts <- read_delim("~/BTSync/FetalRNAseq/Github/GENEX-FB1/Results/MvsF_12_20_PCW_FDR_0.1/tables/MalevsFemale.complete.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  filter(! is.na(padj)) %>%
  dplyr::select(Id, starts_with('norm'))
colnames(counts) <- str_replace(colnames(counts), 'norm\\.', '')
counts <- counts[ , order(names(counts))]
counts <- dplyr::select(counts, id=Id, everything())

# import co-variates
# Seems that data are coerced into float/integer column-wiseduring import with 
# the LoadFile command even though the data are encoded row-wise. I think the 
# only solution to this is to round every thing to whole numbers.
target <- read_delim("~/BTSync/FetalRNAseq/Github/GENEX-FB1/Shiny/GENEX-FB1/Data/target.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(label=as.character(label))
target$Sex<-as.integer(factor(target$Sex))
target$Sequencer<-as.integer(factor(target$Sequencer))
target$RIN<- round(target$RIN, digits=0)
target <- arrange(target, label)
target <- as.data.frame(target)
rownames(target) <- target$label
target <- dplyr::select(target, -files, -label, -`Read length`)
target<- as.data.frame(t(target))
target$id=rownames(target)
target <- target[ , order(names(target))]
target <- dplyr::select(target, id, everything())

# convert vcf file to table of genotypes and list of SNP locations
individual_vcf_files = "/Users/heo3/BTSync/FetalRNAseq/Imputation/chr*.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.gz"
vcf_file_name = "/Users/heo3/BTSync/FetalRNAseq/Imputation/chr_all.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.gz"
python_script_name = "~/BTsync/FetalRNAseq/LabNotes/Python/VCF2numeric.py"
snps_location_file_name = "~/BTSync/FetalRNAseq/Github/GENEX-FB2/MatrixEQTL/snp_pos.txt"
snps_file_name = "~/BTSync/FetalRNAseq/Github/GENEX-FB2/MatrixEQTL/genotypes.txt"

# use file tests to avoid re-doing unnecessary shell commands, which are pretty time consuming
if (! file_test("-f", vcf_file_name)) {
  system(paste('bcftools concat', individual_vcf_files, '>', vcf_file_name))
}

if (! file_test("-f", snps_location_file_name) | ! file_test("-f", snps_file_name)) {
  system(paste('bcftools view', vcf_file_name, '| python', python_script_name, snps_file_name, snps_location_file_name))
}

snp_tbl<- read_delim(snps_file_name, delim='\t')
snp_tbl <- dplyr::rename(snp_tbl,  `18208` = `18121`)
snp_tbl <- snp_tbl[ , order(names(snp_tbl))]
snp_tbl <- dplyr::select(snp_tbl, id, everything())

# remove samples that are missing in one or more of the files
counts <- counts %>% dplyr::select(-one_of(setdiff(colnames(counts), colnames(target))))
counts <- counts %>% dplyr::select(-one_of(setdiff(colnames(counts), colnames(snp_tbl))))
target <- target %>% dplyr::select(-one_of(setdiff(colnames(target), colnames(counts))))
snp_tbl <- snp_tbl %>% dplyr::select(-one_of(setdiff(colnames(snp_tbl), colnames(counts))))

# save modified files and read in as matrix EQTL objects
counts_file_name = tempfile()
write_tsv(counts, counts_file_name)

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(counts_file_name);


covariates_file_name = tempfile()
write_tsv(target, covariates_file_name)

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

geneotype_file_name = tempfile()
write_tsv(snp_tbl, geneotype_file_name)

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( geneotype_file_name );

#extract gene positions from GTF file
gene_location_file_name = "~/BTSync/FetalRNAseq/Github/GENEX-FB2/MatrixEQTL/geneloc.txt"
system(paste('echo "geneid\tchr\ts1\ts2" >', gene_location_file_name))
system(paste("cat ~/BTSync/FetalRNAseq/Ref/genes.gtf | awk '{if ($3 == \"gene\") print $10, $1, $4, $5}' | sed 's/[\";]//g' >>", gene_location_file_name))
gene_location_file_name = "~/BTSync/FetalRNAseq/Github/GENEX-FB2/MatrixEQTL/geneloc.txt"
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);


###############################################################
output_file_name_cis = "~/BTSync/FetalRNAseq/Github/GENEX-FB2/MatrixEQTL/cis_eqtl.txt"
output_file_name_tra = "~/BTSync/FetalRNAseq/Github/GENEX-FB2/MatrixEQTL/trans_eqtl.txt"

pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 1e-2;

cisDist = 1e6;

useModel = modelLINEAR

errorCovariance = numeric()

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

write_tsv(filter(me$cis$eqtls, FDR <=.1), "~/BTSync/FetalRNAseq/Github/GENEX-FB2/MatrixEQTL/cis_eqtl.txt")
write_tsv(filter(me$trans$eqtls, FDR <=.1), "~/BTSync/FetalRNAseq/Github/GENEX-FB2/MatrixEQTL/trans_eqtl.1.txt")
write_tsv(filter(me$trans$eqtls, FDR <=.1e-10), "~/BTSync/FetalRNAseq/Github/GENEX-FB2/MatrixEQTL/trans_eqtl.txt")

save.image(file="~/BTSync/FetalRNAseq/Github/GENEX-FB2/MatrixEQTL/results.RData")