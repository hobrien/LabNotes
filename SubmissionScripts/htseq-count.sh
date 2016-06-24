#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
filename=${1##*/}
sampleID=${filename%%.bam}

samtools sort -n /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.bam \
  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.nsort 
  
htseq-count -f bam -s reverse -t exon -i gene_id -m intersection-strict \
  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.nsort.bam \
  /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
  > /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.counts.txt
