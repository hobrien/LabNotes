#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
outfile1=$(basename $1 .fastq.gz)_trimmed.fastq.gz
outfile2=$(basename $2 .fastq.gz)_trimmed.fastq.gz

/home/mpmho/.local/bin/cutadapt \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
            -o /home/mpmho/Trimmed/$outfile1 -p /home/mpmho/Trimmed/$outfile2 $@

