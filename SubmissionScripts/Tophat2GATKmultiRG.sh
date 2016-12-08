#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

folder=$1
SampleID=$2

echo "extracting reads from $SampleID-1"
if [[ `grep -P "\s$SampleID-2(\s|$)" ~/LabNotes/sequences.txt | wc -l` != 2 ]]
then
    echo "incorrect number of rows for $SampleID-2"
    exit 1
fi
rgid=`find /c8000xd3/databank/foetal-rna/ /c8000xd2/foetalRNAseq/ -name $(grep -P "\s$SampleID-2(\s|$)"  ~/LabNotes/sequences.txt | head -1 | cut -f 1)*fastq.gz | xargs zcat | head -1 | cut -d: -f 3,4'`
if [ ! $rgid ]
then
    echo "trying without gzip"
    rgid=`find /c8000xd3/databank/foetal-rna/ /c8000xd2/foetalRNAseq/ -name $(grep -P "\s$SampleID-2(\s|$)"  ~/LabNotes/sequences.txt | head -1 | cut -f 1)*fastq | xargs head -1 | cut -d: -f 3,4'`
fi
if [ ! $rgid ]
then
    echo "Could not get RGID from $SampleID-2"
    exit 0
fi
~/bin/samtools view -h $folder/accepted_hits_fixup_merge_sort.bam | grep -v $rgid | samtools view -bhS > $folder/accepted_hits_fixup_merge_sort_1.bam
if [ $? -eq 0 ]
then
    echo "Finished extracting reads from $SampleID-1"
else
    echo "could not extract reads from $SampleID-1"
    exit 1
fi

echo "Adding read groups for $SampleID-1"
bash ~/LabNotes/SubmissionScripts/AddRG.sh $folder/accepted_hits_fixup_merge_sort_1.bam $SampleID-1
if [ $? -eq 0 ]
then
    echo "Finished adding read groups for $SampleID-1"
else
    echo "could not add read groups for $SampleID-1"
    exit 1
fi


echo "extracting reads from $SampleID-2"
if [[ `grep -P "\s$SampleID-1(\s|$)" ~/LabNotes/sequences.txt | wc -l` != 2 ]]
then
    echo "incorrect number of rows for $SampleID-1"
    exit 1
fi
rgid=`find /c8000xd3/databank/foetal-rna/ /c8000xd2/foetalRNAseq/ -name $(grep -P "\s$SampleID-1(\s|$)"  ~/LabNotes/sequences.txt | head -1 | cut -f 1)*fastq.gz | xargs zcat | head -1 | cut -d: -f 3,4'`
if [ ! $rgid ]
then
    echo "trying without gzip"
    rgid=`find /c8000xd3/databank/foetal-rna/ /c8000xd2/foetalRNAseq/ -name $(grep -P "\s$SampleID-1(\s|$)"  ~/LabNotes/sequences.txt | head -1 | cut -f 1)*fastq | xargs head -1 | cut -d: -f 3,4'`
fi
if [ ! $rgid ]
then
    echo "Could not get RGID from $SampleID-1"
    exit 0
fi
~/bin/samtools view -h $folder/accepted_hits_fixup_merge_sort.bam | grep -v $rgid | samtools view -bhS > $folder/accepted_hits_fixup_merge_sort_2.bam
if [ $? -eq 0 ]
then
    echo "Finished extracting reads from $SampleID-2"
else
    echo "could not extract reads from $SampleID-2"
    exit 1
fi

echo "Adding read groups for $SampleID-2"
bash ~/LabNotes/SubmissionScripts/AddRG.sh $folder/accepted_hits_fixup_merge_sort_2.bam $SampleID-2
if [ $? -eq 0 ]
then
    echo "Finished adding read groups for $SampleID-2"
else
    echo "could not add read groups for $SampleID-2"
    exit 1
fi

echo "Merging $Sample-1 with $SampleID-2"
bash ~/LabNotes/SubmissionScripts/MergeSAM.sh $folder/accepted_hits_fixup_merge_sort_1_RG.bam $folder/accepted_hits_fixup_merge_sort_2_RG.bam
if [ $? -eq 0 ]
then
    echo "Finished merging $Sample-1 with $SampleID-2"
else
    echo "could not merging $Sample-1 with $SampleID-2"
    exit 1
fi
mv $folder/accepted_hits_fixup_merge_sort_1_RG_merge.bam $folder/accepted_hits_fixup_merge_sort_RG.bam