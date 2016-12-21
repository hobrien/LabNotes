#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

java -jar ~/src/picard-tools-1.139/picard.jar MergeSamFiles \
      I=$1 \
      I=$2 \
      O=${1%.*}_merge.bam 
exit $?

