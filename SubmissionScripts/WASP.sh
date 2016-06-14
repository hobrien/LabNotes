#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
#path=${1%/*}
#sampleID=${path##*/}
python ~/src/WASP/mapping/find_intersecting_snps.py -p $1 /c8000xd3/rnaseq-heath/Genotypes/Imputation2/SNPs/
