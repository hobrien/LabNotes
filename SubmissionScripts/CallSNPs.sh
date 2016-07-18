#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH


samtools mpileup -uf  /home/mpmho/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa $@ \
   | bcftools view -bvcg - > /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/var.raw.bcf  
   
bcftools view /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/var.raw.bcf \
   | vcfutils.pl varFilter -D1000 > /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/var.flt.vcf 
