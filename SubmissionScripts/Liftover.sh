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

~/bin/java -Xmx20g -jar ~/src/picard-tools-2.4.1/picard.jar LiftoverVcf \
  INPUT=/c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr1.recoded100k.vcf.gz \
  OUTPUT=/c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr1.GRCh38.vcf.gz \
  CHAIN=/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/hg19ToHg38.over.chain \
  REJECT=/c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr1.reject.vcf.gz \
  REFERENCE_SEQUENCE=/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa 
  #MAX_RECORDS_IN_RAM=5000000 
