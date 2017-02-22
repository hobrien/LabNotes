#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution

echo "Getting counts for $filename"
  
htseq-count -f bam -s reverse -t exon -i gene_id -m intersection-strict \
  $1 \
  /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
  > ${1%%.bam}_counts.txt
echo "Finished getting counts for $filename"

