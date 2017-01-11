#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

BASEDIR=/c8000xd3/rnaseq-heath/Mappings/
# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
filename=${1##*/}
SampleID=${filename%%_*} #This will remove the library ID, lane number and read number

# deal with duplicated samples
# this will keep incrementing the SampleID until a unique one is found
replicate=1
baseID=$SampleID
SampleID=$baseID-$replicate
while [[ -d $BASEDIR/$SampleID ]]
do
  replicate=$((replicate+1))
  SampleID=$baseID-$replicate
done
  
echo "Starting mapping for $SampleID"
mkdir /c8000xd3/rnaseq-heath/Mappings/$SampleID
if [ ! -f $BASEDIR/$SampleID/BAM/accepted_hits.bam ] || [ ! -f $BASEDIR/$SampleID/BAM/unmapped.bam ]
then
    echo "Tophat mapping for $SampleID"
    tophat --keep-fasta-order --library-type fr-secondstrand --mate-inner-dist 500 --mate-std-dev 50 --num-threads 8 \
      --transcriptome-index /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.inx \
      --output-dir $BASEDIR/$SampleID \
      /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome $@
    if [ $? -eq 0 ]
    then
        echo "Finished Tophat mapping for $SampleID"
    else
        echo "Tophat mapping failed on $SampleID"
        exit 1
    fi    
    mkdir /c8000xd3/rnaseq-heath/Mappings/$SampleID/BAM
    mv $BASEDIR/$SampleID/accepted_hits.bam $BASEDIR/$SampleID/BAM/
    mv $BASEDIR/$SampleID/unmapped.bam $BASEDIR/$SampleID/BAM/
fi

if [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.sort.bam ]
then
    echo "Sorting $SampleID"
    samtools sort $BASEDIR/$SampleID/BAM/accepted_hits.bam $BASEDIR/$SampleID/BAM/$SampleID.sort
    if [ $? -eq 0 ]
    then
        echo "Finished sorting for $SampleID"
    else
        echo "Could not sort BAM for $SampleID"
        exit 1
    fi
fi   

if [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.sort.bam.bai ]
then
    echo "Indexing $SampleID"
    samtools index $BASEDIR/$SampleID/BAM/$SampleID.sort.bam    
    if [ $? -eq 0 ]
    then
        echo "Finished indexing for $SampleID"
    else
        echo "Could not index BAM for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.chr.bam ]
then
    echo "Running RNAseqQC $SampleID"
    bash ~/LabNotes/SubmissionScripts/RNAseqQC.sh $BASEDIR/$SampleID/BAM/$SampleID.sort.bam   
    if [ $? -eq 0 ]
    then
        echo "Finished running RNAseqQC for $SampleID"
    else
        echo "Could not run RNAseqQC for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.chr.counts.txt ]
then
    echo "Running htseq-count $SampleID"
    bash ~/LabNotes/SubmissionScripts/htseq-count.sh $BASEDIR/$SampleID/BAM/$SampleID.chr.bam   
    if [ $? -eq 0 ]
    then
        echo "Finished running htseq-count for $SampleID"
    else
        echo "Could not run htseq-count for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/BAM/accepted_hits_filtered_sort_dedup.bam ]
then
    echo "Running WASP remapping on $SampleID"
    bash ~/LabNotes/SubmissionScripts/WASP.sh $BASEDIR/$SampleID/BAM/accepted_hits.bam
    if [ $? -eq 0 ]
    then
        echo "Finished running WASP remapping for $SampleID"
    else
        echo "Could not run WASP remapping for $SampleID"
        exit 1
    fi
fi   


if [ ! -f $BASEDIR/$SampleID/BAM/accepted_hits_filtered_sort_dedup.bam ] || [ ! -f $BASEDIR/$SampleID/BAM/unmapped_fixup.bam ]
then
    echo "Fixing BAM formatting for $SampleID"
    bash ~/LabNotes/SubmissionScripts/tophat-recondition.sh $BASEDIR/$SampleID/BAM/
    if [ $? -eq 0 ]
    then
        echo "Finished BAM reformatting on $SampleID"
    else
        echo "tophat-recondition failed on $SampleID"
        exit 1
    fi    
fi


if [ ! -f $BASEDIR/$SampleID/BAM/accepted_hits_filtered_sort_dedup_sort.bam ]
then
    echo "Sorting BAM for $SampleID"
    # Sort BAM files by query name
    bash ~/LabNotes/SubmissionScripts/SamtoolsSort.sh $BASEDIR/$SampleID/BAM/accepted_hits_filtered_sort_dedup.bam
    if [ $? -eq 0 ]
    then
        echo "Finished BAM sorting on $SampleID"
    else
        echo "Samtools sort failed on $SampleID"
        exit 1
    fi    
fi

if [ ! -f $BASEDIR/$SampleID/BAM/accepted_hits_filtered_sort_dedup_sort_RG.bam ]
then
    echo "Adding read groups for $SampleID"
    bash ~/LabNotes/SubmissionScripts/AddRG.sh $BASEDIR/$SampleID/BAM/accepted_hits_filtered_sort_dedup_sort.bam $SampleID
    if [ $? -eq 0 ]
    then
        echo "Finished adding read groups for $SampleID"
    else
        echo "could not add read groups for $SampleID"
        exit 1
    fi
fi


#bash ~/LabNotes/SubmissionScripts/dexseq-count.sh /c8000xd3/rnaseq-heath/Mappings/$SampleID/BAM/$SampleID.chr.bam
#bash ~/LabNotes/SubmissionScripts/DivideBAM.sh $SampleID
#bash ~/LabNotes/SubmissionScripts/CallSNPs.sh /c8000xd3/rnaseq-heath/Mappings/$SampleID/BAM/Chromosomes/$SampleID.chr22.bam
#bash ~/LabNotes/SubmissionScripts/GTcheck.sh $SampleID
#index=`grep $SampleID ~/LabNotes/VCFindex.txt | cut -f 2`
#bash ~/LabNotes/SubmissionScripts/WASPnonRef.sh $SampleID $index
#bash ~/LabNotes/SubmissionScripts/RNAseqQCwasp.sh /c8000xd3/rnaseq-heath/Mappings/$SampleID/BAM/$SampleID.chr.nonref.merged.sorted.bam
#samtools sort /c8000xd3/rnaseq-heath/Mappings/$SampleID/BAM/$SampleID.chr.nonref.merged.dedup.bam /c8000xd3/rnaseq-heath/Mappings/$SampleID/BAM/$SampleID.chr.nonref.merged.dedup.sort
#bash ~/LabNotes/SubmissionScripts/RNAseqQCwasp.sh /c8000xd3/rnaseq-heath/Mappings/$SampleID/BAM/$SampleID.chr.nonref.merged.dedup.sort.bam
#bash ~/LabNotes/SubmissionScripts/clipOverlap.sh /c8000xd3/rnaseq-heath/Mappings/$SampleID/BAM/$SampleID.chr.nonref.merged.dedup.sort.bam
#samtools index /c8000xd3/rnaseq-heath/Mappings/$SampleID/BAM/$SampleID.chr.nonref.merged.dedup.sort.clip.bam
