#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
echo "running fastQC on $@"
function joinStrings { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }
#zcat '@$' | ~/src/FastQC/fastqc --outdir=/home/mpmho/FastQC stdin
~/src/FastQC/fastqc --outdir=/home/mpmho/FastQC "$@"

echo "finished running fastQC on $@"
