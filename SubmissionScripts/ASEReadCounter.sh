#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
filename=${1##*/}
sampleID=${filename%%_*} #This will remove the library ID, lane number and read number
chr=$2

java -jar /share/apps/GenomeAnalysisTK.jar -T ASEReadCounter \
   -R /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa \
   -o ${1%.*}.$chr.ase.rtable \
   -I $1 \
   -sites /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.filter_dup.vcf.gz \
   -U ALLOW_N_CIGAR_READS \
   -minDepth 10 \
   --minMappingQuality 10 \
   --minBaseQuality 2
exit $?
