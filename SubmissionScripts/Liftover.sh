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
bcftools view /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$1.dose.vcf.gz \
  |perl -pe 's/^(\d+)/chr$1/' \
  | bgzip \
  > /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$1.recoded.vcf.gz
CrossMap.py vcf /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/hg19ToHg38.over.chain \
  /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$1.recoded.vcf.gz \
  /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa \
  /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$1.GRCh38.vcf
bgzip /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$1.GRCh38.vcf
rm /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$1.GRCh38.vcf 
rm /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$1.recoded.vcf.gz
