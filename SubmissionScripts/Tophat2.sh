#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
filename=${1##*/}
sampleID=${filename%%_L00*} #This will remove the lane number and read number
echo "Starting mapping for $sampleID"
mkdir /c8000xd3/rnaseq-heath/Mappings/$sampleID
tophat --keep-fasta-order --library-type fr-secondstrand --mate-inner-dist 500 --mate-std-dev 50 --num-threads 8 \
  --transcriptome-index /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.inx \
  --output-dir /c8000xd3/rnaseq-heath/Mappings/$sampleID \
  /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome $@
mkdir /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM
mv /c8000xd3/rnaseq-heath/Mappings/$sampleID/accepted_hits.bam /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/
mv /c8000xd3/rnaseq-heath/Mappings/$sampleID/unmapped.bam /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/
echo "Sorting and indexing $sampleID"
samtools sort /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/accepted_hits.bam /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.sort
samtools index /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.sort.bam
echo "Finished mapping for $sampleID"
bash ~/LabNotes/SubmissionScripts/RNAseqQC.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.sort.bam
bash ~/LabNotes/SubmissionScripts/htseq-count.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.bam
bash ~/LabNotes/SubmissionScripts/dexseq-count.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.bam
bash ~/LabNotes/SubmissionScripts/DivideBAM.sh $sampleID
bash ~/LabNotes/SubmissionScripts/CallSNPs.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/Chromosomes/$sampleID.sort.bam
bash ~/SubmissionScripts/GTcheck.sh $sampleID

