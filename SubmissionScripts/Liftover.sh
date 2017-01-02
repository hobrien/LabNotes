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
#bcftools view /c8000xd3/rnaseq-heath/Genotypes/Imputation3/hg19/chr$chr.dose.rename.filter_samples.filter_sites.rsID.vcf.gz \
#  | perl -pe 's/^(##contig=<ID=)?(\d+)/$1chr$2/' \
#  | bgzip \
#  > /c8000xd3/rnaseq-heath/Genotypes/Imputation3/hg19/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.vcf.gz
#CrossMap.py vcf /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/hg19ToHg38.over.chain \
#  /c8000xd3/rnaseq-heath/Genotypes/Imputation3/hg19/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.vcf.gz \
#  /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Alt/genome.fa \
#  /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.vcf
#bgzip /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.vcf
#bash ~/LabNotes/SubmissionScripts/SortVCF.sh /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.vcf.gz
#bash ~/LabNotes/SubmissionScripts/checkVCF.sh /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.gz

cut -f 2 /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.check.nonSnp > /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.excludedSNPs.txt
cut -f 2 /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.check.ref | cut -f 2 -d: >> /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.excludedSNPs.txt
exset=`sort /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.excludedSNPs.txt |uniq | perl -pe 's/^/POS=/' | paste -s -d'|'`
if [[ ${#exset} > 0 ]]
then
    bcftools filter -Oz -e $exset /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.gz \
      > /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter.vcf.gz
    bash ~/LabNotes/SubmissionScripts/checkVCF.sh /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter.vcf.gz
else
    cp /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.gz \
      /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter.vcf.gz
fi
bcftools index /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter.vcf.gz
exit $?
