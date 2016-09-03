#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

refID=${1%-2}
echo $refID > /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/$1.txt
bcftools reheader -s /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/$1.txt -o /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/$1.chr22.renamed.bcf /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/$1.chr22.raw.bcf 
rm /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/$1.txt 
bcftools index /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/$1.chr22.renamed.bcf
bcftools gtcheck -g /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr22.GRCh38.renamed.vcf.gz -p chr22_$1 -s $1 -S $refID /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/$1.chr22.renamed.bcf
