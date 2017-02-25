#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=20G
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

infile=$1
outfile=${1%.*}_RG.bam
SampleID=$2
if [[ `grep -P "\s$SampleID(\s|$)" ~/LabNotes/sequences.txt | wc -l` != 2 ]]
then
    echo "incorrect number of rows for $SampleID"
    exit 0
fi
rgid=`find /c8000xd2/foetalRNAseq/ /c8000xd3/databank/foetal-rna/ -name $(grep -P "\s$SampleID(\s|$)"  ~/LabNotes/sequences.txt | head -1 | cut -f 1)*fastq.gz | xargs zcat | head -1 | cut -d: -f 3,4,10 | perl -pe 's/:/./g'`
if [ ! $rgid ]
then
    echo "trying without gzip"
    rgid=`find /c8000xd2/foetalRNAseq/ /c8000xd3/databank/foetal-rna/ -name $(grep -P "\s$SampleID(\s|$)"  ~/LabNotes/sequences.txt | head -1 | cut -f 1)*fastq | xargs head -1 | cut -d: -f 3,4,10 | perl -pe 's/:/./g'`
fi
if [ ! $rgid ]
then
    echo "Could not get RGID from $SampleID"
    exit 0
fi

java -Xmx2g -jar ~/src/picard-tools-1.139/picard.jar AddOrReplaceReadGroups \
      I=$infile \
      O=$outfile \
      RGID=$rgid \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=$rgid \
      RGSM=$SampleID
exit $?
