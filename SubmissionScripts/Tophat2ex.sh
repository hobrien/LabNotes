#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
filename1=${1##*/}
filename2=${2##*/}

if [ $filename1 == 15533_2_TGACCA_L008_R1_001.fastq ]
then
    sampleID=15533_2
else    
    sampleID=${filename1%%_*}
fi

if [ $sampleID == 16810 ] || [ $sampleID == 15533 ]
then
    read1=/c8000xd3/rnaseq-heath/Trimmed/${filename1%%.*}_val_1.fq.gz
    read2=/c8000xd3/rnaseq-heath/Trimmed/${filename2%%.*}_val_2.fq.gz
else
    read1=/c8000xd3/rnaseq-heath/Trimmed/${filename1%%.*}_val_1.fq
    read2=/c8000xd3/rnaseq-heath/Trimmed/${filename2%%.*}_val_2.fq
fi
echo "Started processing $sampleID"
echo $read1
echo $read2
mkdir /c8000xd3/rnaseq-heath/Mappings/$sampleID
tophat --keep-fasta-order --library-type fr-secondstrand --mate-inner-dist 500 --mate-std-dev 50 --num-threads 8 \
  --transcriptome-index /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.inx \
  --output-dir /c8000xd3/rnaseq-heath/Mappings/$sampleID \
  /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome \
  $read1 $read2
mkdir /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM
mv /c8000xd3/rnaseq-heath/Mappings/$sampleID/accepted_hits.bam /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/
mv /c8000xd3/rnaseq-heath/Mappings/$sampleID/unmapped.bam /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/
samtools sort /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/accepted_hits.bam /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.sort
samtools index /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.sort.bam
echo "Finished processing $sampleID"
