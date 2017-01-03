#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

# Usage: qsub SortSam.sh FILENAME SORT_ORDER
# possible values for SORT_ORDER: unsorted, queryname, coordinate, duplicate

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH



samtools merge ${1%.*}_merge.bam $@
exit $?
