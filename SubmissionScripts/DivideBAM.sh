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
    folder=${dataset##*/}
    chr=22
    echo $folder_path/BAM/$folder.$chr.bam
    #samtools view -bh $dataset $chr > $folder_path/BAM/$folder.chr.bam
    #samtools index /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.bam
done
