#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# Index is the position of the sample in the VCF plus 9 (ie; sample 1 has an index of 10) 
# Name is the BrainBankID
sampleName=$1
sampleIndex=$2

mkdir  /c8000xd3/rnaseq-heath/Genotypes/Imputation2/$sampleName/
for num in {1..21}
do
    bcftools view -H  /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$num.GRCh38.vcf.gz | cut -f1,2,3,4,5,6,7,8,9,$sampleIndex | grep -v '0|0' | cut -f2,4,5 |gzip -c > /c8000xd3/rnaseq-heath/Genotypes/Imputation2/$sampleName/chr$num.snps.txt.gz
done
