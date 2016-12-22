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
chr=$1
bcftools view /c8000xd3/rnaseq-heath/Genotypes/Imputation3/hg19/chr$chr.dose.rename.filter_samples.filter_sites.rsID.vcf.gz \
  |perl -pe 's/^(\d+)/chr$1/' \
  | bgzip \
  > /c8000xd3/rnaseq-heath/Genotypes/Imputation3/hg19/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.vcf.gz
CrossMap.py vcf /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/hg19ToHg38.over.chain \
  /c8000xd3/rnaseq-heath/Genotypes/Imputation3/hg19/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.vcf.gz \
  /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa \
  /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.vcf
bgzip /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.vcf
rm /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.vcf
exit $?