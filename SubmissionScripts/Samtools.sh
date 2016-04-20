#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

folder_path=${@%/*}
folder=${folder_path##*/}
#samtools sort -o $folder_path/$folder.sort.bam $@
samtools index $folder_path/$folder.sort.bam
#samtools mpileup -d 8000 -f /home/mpmho/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa -r chr1:98046571-98046571 $@
