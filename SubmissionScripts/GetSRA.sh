#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
#path=${1%/*}
sraID=$1

echo "Starting SRA dump for $sraID"
cd /c8000xd3/rnaseq-heath/SRA
fastq-dump --split-files --gzip --readids $sraID
echo "Finished SRA dump for $sraID. Output in /c8000xd3/rnaseq-heath/SRA"
