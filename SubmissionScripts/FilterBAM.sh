#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
sampleID=$1

#mkdir /c8000xd3/rnaseq-heath/Cufflinks/$sampleID

#samtools sort /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.ex.bam /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.ex.sort
samtools index /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.ex.bam
samtools view -bh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.ex.bam \
chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
> /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.bam
#rm /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.ex.sort.bam /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.ex.sort.bam.bai
#samtools view -H /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.bam |grep -E '@HD|@PG|chr[0-9]{1,2}\s' > /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/header.sam
#samtools reheader /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/header.sam /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.bam \
#> /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.reheader.bam
samtools index /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.bam


#cufflinks -p 8 -o /c8000xd3/rnaseq-heath/Cufflinks/$sampleID --library-type fr-secondstrand \
#  -g /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
#  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.bam
