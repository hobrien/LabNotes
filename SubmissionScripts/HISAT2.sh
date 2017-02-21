#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 8
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH
REF=/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta
ANNOTATION=/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode
MAPPER=HISAT2
# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
SampleID=$1 
BASEDIR=/c8000xd3/rnaseq-heath/Mappings/$SampleID

echo "Starting mapping for $SampleID"

if [ ! -d $BASEDIR ]
then
    echo "Creating directory $SampleID"
    mkdir $BASEDIR
    if [ $? -eq 0 ]
    then
        echo "Finished creating directory for $SampleID"
    else
        echo "Could not create directory for $SampleID"
        exit 1
    fi    
fi


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

if [ ! -f $BASEDIR/BAM/accepted_hits.bam ]
then
    echo "$MAPPER mapping for $SampleID"
    sequences=$(for name in `grep -P "\s$SampleID(\s|$)"  ~/LabNotes/sequences.txt | cut -f 1`; do find /c8000xd2/foetalRNAseq/ /c8000xd3/databank/foetal-rna/ /c8000xd3/rnaseq-heath/SRA -name $name*fastq*; done)
    hisat2 --fr --threads 8 -x $REF/genome  --known-splicesite-infile $ANNOTATION/splicesites.txt \
      -1 $sequences \
      | samtools view -S -bo $BASEDIR/accepted_hits.bam -
    if [ $? -eq 0 ]
    then
        echo "Finished $MAPPER mapping for $SampleID"
    else
        echo "$MAPPER mapping failed for $SampleID"
        exit 1
    fi    
    mkdir $BASEDIR/BAM
    mv $BASEDIR/accepted_hits.bam $BASEDIR/BAM/
fi

if [ ! -f $BASEDIR/BAM/$SampleID.sort.bam ]
then
    echo "Sorting $SampleID"
    samtools sort $BASEDIR/BAM/accepted_hits.bam $BASEDIR/BAM/$SampleID.sort
    if [ $? -eq 0 ]
    then
        echo "Finished sorting for $SampleID"
    else
        echo "Could not sort BAM for $SampleID"
        exit 1
    fi
fi   

if [ ! -f $BASEDIR/BAM/$SampleID.sort.bam.bai ]
then
    echo "Indexing $SampleID"
    samtools index $BASEDIR/BAM/$SampleID.sort.bam    
    if [ $? -eq 0 ]
    then
        echo "Finished indexing for $SampleID"
    else
        echo "Could not index BAM for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/${SampleID}_stats.txt ]   
then
    echo "Running RSeQC stats on $MAPPER mapping for $SampleID"
   bam_stat.py -i $BASEDIR/BAM/$SampleID.sort.bam > $BASEDIR/${SampleID}_stats.txt
    if [ $? -eq 0 ]
    then
        echo "Finished running RSeQC stats on $MAPPER mapping for $SampleID"
    else
        echo "RSeQC stats could not be run on $MAPPER mapping for $SampleID"
        exit 1
    fi    
fi


if [ ! -f $BASEDIR/BAM/${SampleID}_nsort.bam ]   
then
    echo "Sorting $MAPPER mapping for $SampleID"
    samtools sort -n $BASEDIR/BAM/accepted_hits.bam \
      $BASEDIR/BAM/${SampleID}_nsort 
    if [ $? -eq 0 ]
    then
        echo "Finished sorting $MAPPER mapping for $SampleID"
    else
        echo "Could not sort $MAPPER mapping for $SampleID"
        exit 1
    fi    
fi

if [ ! -f $BASEDIR/${SampleID}_counts.txt ]   
then
    echo "Running htseq-count on sorted $MAPPER mapping for $SampleID"
    htseq-count -f bam -s reverse -t exon -i gene_id -m intersection-strict \
      $BASEDIR/BAM/${SampleID}_nsort.bam \
      /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
      > $BASEDIR/${SampleID}_counts.txt
    if [ $? -eq 0 ]
    then
        echo "Finished running htseq-count on sorted $MAPPER mapping for $SampleID"
    else
        echo "htseq-count could not be run on sorted $MAPPER mapping for $SampleID"
        exit 1
    fi    
fi

