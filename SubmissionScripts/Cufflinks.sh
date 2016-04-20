#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
filename1=${1##*/}
sampleID=${filename1%%.*}

mkdir /c8000xd3/rnaseq-heath/Cufflinks/$sampleID

cufflinks -p 8 -o /c8000xd3/rnaseq-heath/Cufflinks/$sampleID --library-type fr-secondstrand \
  -g /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf $@