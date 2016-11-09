#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
for input in `find /c8000xd3/rnaseq-heath/Mappings/ -name *.chr.nonref.merged.dedup.sort.clip.bam`
do
    filename1=${input##*/}
    file_location=${input%/*} 
    sampleID=${filename1%%.*}
    echo "Started processing $sampleID"



    samtools index $input
    echo "Finished processing $sampleID"
done
