#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

java -Xmx2g -jar ~/src/picard-tools-1.139/picard.jar MergeBamAlignment \
       ALIGNED=$1 \ 
       UNMAPPED=$2 \ 
       O=${1%.*}.withUnmapped.bam \
       R=/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa