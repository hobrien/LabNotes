#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

# Usage: qsub SortSam.sh FILENAME SORT_ORDER
# possible values for SORT_ORDER: unsorted, queryname, coordinate, duplicate

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

java -jar ~/src/picard-tools-1.139/picard.jar SortSam \
  I=$1 \
  O=${1%.*}_$2.bam \
  SORT_ORDER=$2
