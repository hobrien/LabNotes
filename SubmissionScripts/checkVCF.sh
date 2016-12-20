#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#


export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

python ~/src/checkVCF-20140116/checkVCF.py \
  #-r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa \
  -r ~/src/checkVCF-20140116/hs37d5.fa \
  -o ${1%.*} $1
