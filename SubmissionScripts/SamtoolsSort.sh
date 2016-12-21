#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

~/bin/samtools sort -o ${1%.*}_sort.bam $@
exit $?
