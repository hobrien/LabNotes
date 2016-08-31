#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution

sampleID=$1
sampleIndex=$2

bash ~/LabNotes/SubmissionScripts/ExtractSNPs.sh $1 $2

echo "Starting WASP non-ref Remapping on $sampleID"

mkdir /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/RemapNonRef
cp /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.bam /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/RemapNonRef/$sampleID.chr.bam

python ~/src/WASP/mapping/find_intersecting_snps.py -p /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/RemapNonRef/$sampleID.chr.bam /c8000xd3/rnaseq-heath/Genotypes/Imputation2/$sampleID/

tophat --keep-fasta-order --library-type fr-secondstrand --mate-inner-dist 500  --mate-std-dev 50 --num-threads 8 \
  --transcriptome-index /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.inx \
  --output-dir /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/RemapNonRef \
  /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome \
  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/RemapNonRef/$sampleID.chr.remap.fq1.gz \
  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/RemapNonRef/$sampleID.chr.remap.fq2.gz

python ~/src/WASP/mapping/filter_remapped_reads.py -p \
  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/RemapNonRef/$sampleID.chr.to.remap.bam \
  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/RemapNonRef/accepted_hits.bam \
  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/RemapNonRef/$sampleID.chr.remap.keep.bam \
  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/RemapNonRef/$sampleID.chr.to.remap.num.gz

samtools merge /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.bam \
  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/RemapNonRef/$sampleID.chr.keep.bam \
  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/RemapNonRef/$sampleID.chr.remap.keep.bam
  
samtools sort /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.bam \
  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.sorted
  
samtools index /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.sorted.bam

echo "Finished WASP Non-ref Remapping on $sampleID"

#bash ~/LabNotes/SubmissionScripts/htseq-count.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.bam

python ~/src/WASP/mapping/rmdup_pe.py /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.sorted.bam \
  /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.sorted.dedup.bam

