#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
function joinStrings { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }
#zcat '@$' | ~/src/FastQC/fastqc --outdir=/home/heath/FastQC stdin
~/src/FastQC/fastqc --outdir=/home/heath/FastQC "$@"
