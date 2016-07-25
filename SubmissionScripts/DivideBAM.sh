#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

for dataset in $@
do
    folder_path=${dataset%/*}
    folder=${folder_path%/*}
    folder=${folder##*/}
    chr=22
    mkdir $folder_path/Chromosomes
    #echo $folder_path/Chromosomes/$folder.chr$chr.bam
    samtools view -bh $dataset chr$chr > $folder_path/Chromosomes/$folder.chr$chr.bam
    samtools index $folder_path/Chromosomes/$folder.chr$chr.bam
done
