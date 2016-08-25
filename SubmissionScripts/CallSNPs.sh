#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# this command will call SNPs on all smaple (except 15533 dup):
# find /c8000xd3/rnaseq-heath/Mappings/ -name *chr.nonref.merged.sorted.bam | grep -v SRR | grep -v 50bp | grep -v _ | xargs qsub ../LabNotes/SubmissionScripts/CallSNPs.sh

# Based on: https://github.com/samtools/bcftools/wiki/HOWTOs#mpileup-calling

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
filename=${1##*/}
filename=${filename#*.}
filename=${filename%%.bam}
#~/bin/samtools mpileup -uf /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa /c8000xd3/rnaseq-heath/Mappings/15240/BAM/Chromosomes/15240.chr22.bam | ~/bin/bcftools call -mv -Ob > /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/chr22.raw.bcf
~/bin/samtools mpileup -uf /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa $@ \
   | ~/bin/bcftools call -mv -Ob > /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/$filename.raw.bcf
   
#bcftools view /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/$filename.raw.bcf \
#   | vcfutils.pl varFilter -D1000 > /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/$filename.filt.vcf
