#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
tophat --keep-fasta-order --transcriptome-index /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.inx --library-type fr-secondstrand --mate-inner-dist 150  --mate-std-dev 50 --num-threads 8 --output-dir /home/heath/Mappings/15533_150 /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome /home/heath/Trimmed/15533_TGACCA_L007_R1_001_trimmed.fastq.gz /home/heath/Trimmed/15533_TGACCA_L007_R2_001_trimmed.fastq.gz
