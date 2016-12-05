#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

# This will attempt to string together all of the steps need to produce a BAM thats suitable for analyses with GATK from the output of Tophat

#folder with unmapped.bam and accepted_hits.bam 
folder=$1
sampleID=$2

# Need to fix a bunch of problems with flags in unmapped.bam and change the program ID so that it is different from the one in accepted_hits.bam
# I use a script from https://github.com/cbrueffer/tophat-recondition that I modified to deal with the program ID problem
bash ~/LabNotes/SubmissionScripts/tophat-recondition.sh $folder

# Merge BAM file
bash ~/LabNotes/SubmissionScripts/MergeBAM.sh $folder/accepted_hits_fixup.bam $folder/unmapped_fixup.bam 

# Add Readgroup Info to bam file
bash ~/LabNotes/SubmissionScripts/AddRG.sh $folder/accepted_hits_fixup_merge.bam $sampleID

# Sort BAM files by query name
bash ~/LabNotes/SubmissionScripts/SortSam.sh $folder/accepted_hits_fixup_merge_RG.bam coordinate
