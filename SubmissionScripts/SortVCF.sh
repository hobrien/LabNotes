#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH


for num in {3..22}
do
    /share/apps/vcftools-0.1.14/bin/vcf-sort /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$num.GRCh38.vcf.gz | bgzip > /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$num.GRCh38.sort.vcf.gz
    bcftools view -H /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$num.GRCh38.sort.vcf.gz | grep -e '|1\|1|' | gzip > /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$num.GRCh38.sort.nonref.vcf.gz
done
