#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

/share/apps/vcftools-0.1.14/bin/vcf-sort $1 | bgzip > ${1%.*}.sort.vcf.gz
