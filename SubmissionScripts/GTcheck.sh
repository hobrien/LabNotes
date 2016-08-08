#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

bcftools gtcheck -g /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr22.GRCh38.renamed.vcf.gz -p chr22_$1 -s $1 -S $1 /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/chr22.renamed.bcf