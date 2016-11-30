#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

~/src/WASP-0.2.1/snp2h5/snp2h5 \
  -- chrom /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/chromInfo.txt
  --format vcf \
  --haplotype haplotypes.h5 \
  --snp_index snp_index.h5 \
  --snp_tab   snp_tab.h5 \
  chr*.GRCh38.sort.vcf.gz