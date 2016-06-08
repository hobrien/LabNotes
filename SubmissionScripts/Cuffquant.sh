#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
#path=${1%/*}
sampleID=$1

echo "Starting Cuffquant on $sampleID"

cuffquant -p 8 -u -library-type fr-secondstrand \
  -o /c8000xd3/rnaseq-heath/Cufflinks/Cuffquant \
  -b /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa \
  /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
  /c8000xd3/rnaseq-heath/Cufflinks/Cuffmerge/merged.gtf \
  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.bam
echo "Finished Cuffquant on $sampleID"
