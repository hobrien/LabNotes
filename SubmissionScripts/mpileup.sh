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



    samtools mpileup -d 8000 -r $1 \
    -f /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa \
    $input |python ~/LabNotes/Python/CountBases.py $sampleID >> /c8000xd3/rnaseq-heath/Counts/$1.clip.counts.txt
    echo "Finished processing $sampleID"
done
