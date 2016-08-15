#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution

for sampleID in `ls /c8000xd3/rnaseq-heath/Mappings`
do
    if [[ $sampleID != "15236"* && $sampleID != "SRR"* ]]
    then
        echo $sampleID
        qsub ../LabNotes/SubmissionScripts/dexseq-count.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nsort.bam
    fi
done
