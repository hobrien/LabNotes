#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

folder_path=${@%/*}
folder=${folder_path##*/}
filename=${1%%.bam}
~/bin/samtools sort -o $filename_sort.bam $@
