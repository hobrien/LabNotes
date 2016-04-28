#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
path=${1%/*}
sampleID=${path##*/}

echo "Starting Cuffcompare on $sampleID"
cuffcompare -V -o /c8000xd3/rnaseq-heath/Cufflinks/$sampleID/$sampleID \
  -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf $@
echo "Finished Cuffcompare on $sampleID"
