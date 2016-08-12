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
filename=${filename%%.bam}
sampleID=${filename%%.*}

echo "Getting counts for $filename"

dexseq_count.py -p yes -r name -s reverse -f bam \
  $1 \
  /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gff \
  > /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$filename.dex_counts.txt
echo "Finished getting counts for $filename"