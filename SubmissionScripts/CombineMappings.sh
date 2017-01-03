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
#mkdir -p $BASEDIR/$SampleID/BAM
#Rscript ~/LabNotes/R/CombineCounts.R $SampleID

bash ~/LabNotes/SubmissionScripts/SamtoolsMerge.sh `ls $BASEDIR | grep $SampleID-`
