#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
filename1=${1##*/}
file_location = ${1%/*} 
sampleID=${filename1%%.*}



samtools mpileup -d 8000 \
-f /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa \
$@ |python ~/LabNotes/Python/CountBases.py > $file_location/$sampleID.mpileup.txt
