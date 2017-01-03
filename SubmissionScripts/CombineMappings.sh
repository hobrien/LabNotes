#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
SampleID=$1
#bash ~/LabNotes/SubmissionScripts/RenameFiles $SampleID
mkdir -p /c8000xd3/rnaseq-heath/Mappings/$SampleID/BAM
Rscript ~/LabNotes/R/CombineCounts.R $SampleID
