#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# this command will call SNPs on all smaple (except 15533 dup):
# find /c8000xd3/rnaseq-heath/Mappings/ -name *chr.nonref.merged.sorted.bam | grep -v SRR | grep -v 50bp | grep -v _ | xargs qsub ../LabNotes/SubmissionScripts/CallSNPs.sh
~/bin/samtools mpileup -uf /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa $@ \
   | bcftools call -mv -o b - > /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/var.raw.bcf  
   
bcftools view /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/var.raw.bcf \
   | vcfutils.pl varFilter -D1000 > /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/var.flt.vcf 
