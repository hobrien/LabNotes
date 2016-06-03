#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
#path=${1%/*}
#sampleID=${path##*/}
bcftools view -H  /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$1.GRCh38.vcf.gz | cut -f2,4,5 |gzip -c > /c8000xd3/rnaseq-heath/Genotypes/Imputation2/SNPs/chr$1.snps.txt.gz
