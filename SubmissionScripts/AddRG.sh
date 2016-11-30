#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

infile=$1
outfile=${1%.*}_RG.bam
SampleID=$2
rgid=`find /c8000xd3/databank/foetal-rna/ -name $(grep $SampleID ~/LabNotes/sequences.txt | head -1 | cut -f 1)* | xargs zcat | head -1 | cut -d: -f 3,4,10 | perl -pe 's/:/./g'`

java -Xmx2g -jar ~/src/picard-tools-1.139/picard.jar AddOrReplaceReadGroups \
      I=$infile \
      O=$outfile \
      RGID=$rgid \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=$rgid \
      RGSM=$SampleID
