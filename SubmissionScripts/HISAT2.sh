#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 8
#

REF=/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta
ANNOTATION=/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode
MAPPER=HISAT2



if [ ! -f $REF/genome.1.ht2 ] && [ ! -f $REF/genome.1.ht21 ]
then
    echo "Building index for $MAPPER"
    hisat2-build -p 8 $REF/genome.fa $REF/genome
    if [ $? -eq 0 ]
    then
        echo "Finished building index for $MAPPER"
    else
        echo "Building index for $MAPPER failed"
        exit 1
    fi    
fi

if [ ! -f $ANNOTATION/splicesites.txt ]
then
    echo "Making list of known splice sites"
    hisat2_extract_splice_sites.py $ANNOTATION/genes.gtf > $ANNOTATION/splicesites.txt
    if [ $? -eq 0 ]
    then
        echo "Finished making list of known splice sites"
    else
        echo "Making list of known splice sites failed"
        exit 1
    fi    
fi

hisat2 --fr --threads 8 -x $REF/genome --known-splicesite-infile $ANNOTATION/splicesites.txt \
      -1 $1 -2 $2 | samtools view -S -bo $3
