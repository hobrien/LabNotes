#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

filename=${dataset##*/}
for dataset in $@
do
    file_location=${dataset%/*} 
    filename=${dataset##*/}
    filename=${filename%.*}
    outfile=$file_location/$filename.clip.bam
    #echo $outfile
    bam clipOverlap --stats --in $dataset --out $outfile
done
