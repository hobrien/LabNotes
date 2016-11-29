#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

java -Xmx2g -jar ~/src/picard-tools-1.139/picard.jar AddOrReplaceReadGroups \
      I=$1 \
      O=O=${1%.*}_RG.bam \
      RGID=$3 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=$3 \
      RGSM=$2
