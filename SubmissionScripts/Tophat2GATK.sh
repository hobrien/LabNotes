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
echo "Fixing BAM formatting for $sampleID"
bash ~/LabNotes/SubmissionScripts/tophat-recondition.sh $folder
if [ $? -eq 0 ]
then
    echo "Finished BAM reformatting on $sampleID"
else
    echo "tophat-recondition failed on $sampleID"
    exit 1
fi
# Merge BAM file
echo "Merging mapped and unmapped BAM files for $sampleID"
rm $folder/accepted_hits_fixup_merge.bam
bash ~/LabNotes/SubmissionScripts/SamtoolsMerge.sh $folder/accepted_hits_fixup.bam $folder/unmapped_fixup.bam 
if [ $? -eq 0 ]
then
    echo "Finished BAM merging on $sampleID"
else
    echo "Samtools merge failed on $sampleID"
    exit 1
fi

echo "Sorting BAM for $sampleID"
# Sort BAM files by query name
bash ~/LabNotes/SubmissionScripts/SamtoolsSort.sh $folder/accepted_hits_fixup_merge.bam
if [ $? -eq 0 ]
then
    echo "Finished BAM sorting on $sampleID"
else
    echo "Samtools sort failed on $sampleID"
    exit 1
fi

echo "Adding read group info to merged BAM for $sampleID"
# Add Readgroup Info to bam file
bash ~/LabNotes/SubmissionScripts/AddRG.sh $folder/accepted_hits_fixup_merge_sort.bam $sampleID
if [ $? -eq 0 ]
then
    echo "Finished adding read group info to $sampleID"
else
    echo "AddRG.sh failed on $sampleID"
    exit 1
fi
exit 0

