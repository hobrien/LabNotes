#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

sampleID=$1
bash ~/LabNotes/SubmissionScripts/CallSNPs.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/Chromosomes/$sampleID.chr22.bam
bash ~/LabNotes/SubmissionScripts/GTcheck.sh $sampleID
index=`grep $sampleID ~/LabNotes/VCFindex.txt | cut -f 2`
bash ~/LabNotes/SubmissionScripts/WASPnonRef.sh $sampleID $index
bash ~/LabNotes/SubmissionScripts/RNAseqQCwasp.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.sorted.bam
samtools sort /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.dedup.bam /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.dedup.sort
bash ~/LabNotes/SubmissionScripts/RNAseqQCwasp.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.dedup.sort.bam
bash ~/LabNotes/SubmissionScripts/clipOverlap.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.dedup.sort.bam
