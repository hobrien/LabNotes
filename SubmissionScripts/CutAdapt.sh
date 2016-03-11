#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
/home/heath/.local/bin/cutadapt \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
            -o /home/heath/Trimmed/150429_D00200_0258_BC6UARANXX_4_IL-TP-019_1.trim.fastq.gz \
            -p /home/heath/Trimmed/150429_D00200_0258_BC6UARANXX_4_IL-TP-019_2.trim.fastq.gz \
            /home/heath/Raw/150429_D00200_0258_BC6UARANXX_4_IL-TP-019_1.fastq.gz \
            /home/heath/Raw/150429_D00200_0258_BC6UARANXX_4_IL-TP-019_2.fastq.gz

