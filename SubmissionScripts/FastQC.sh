#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
echo "running fastQC on $@"
#zcat '@$' | ~/src/FastQC/fastqc --outdir=/home/mpmho/FastQC stdin
~/src/FastQC/fastqc --outdir=/home/mpmho/FastQC $@
echo "finished running fastQC on $@"
