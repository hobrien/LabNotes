#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

BASEDIR=/c8000xd3/rnaseq-heath/Mappings/

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
SampleID=$1
#bash ~/LabNotes/SubmissionScripts/RenameFiles $SampleID
mkdir -p $BASEDIR/$SampleID/BAM
Rscript ~/LabNotes/R/CombineCounts.R $SampleID

bash ~/LabNotes/SubmissionScripts/MergeSAM.sh `find $BASEDIR -name accepted_hits_filtered_sort_dedup_sort_RG.bam | grep $SampleID- | sort`
mv $BASEDIR/$SampleID-1/BAM/accepted_hits_filtered_sort_dedup_sort_RG_merge.bam $BASEDIR/$SampleID/BAM/$SampleID.bam
