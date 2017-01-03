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
sampleID=${filename%%_*} #This will remove the library ID, lane number and read number

# deal with duplicated samples
# this will keep incrementing the sampleID until a unique one is found
replicate=1
baseID=$sampleID
while [[ -d $BASEDIR/$sampleID ]]
  do
  replicate=$((replicate+1))
  sampleID=$baseID-$replicate
  done
  
echo "Starting mapping for $sampleID"
mkdir /c8000xd3/rnaseq-heath/Mappings/$sampleID
if [ ! -f $BASEDIR/$sampleID/BAM/accepted_hits.bam ] || [ ! -f $BASEDIR/$sampleID/BAM/unmapped.bam ]
then
    echo "Tophat mapping for $SampleID"
    tophat --keep-fasta-order --library-type fr-secondstrand --mate-inner-dist 500 --mate-std-dev 50 --num-threads 8 \
      --transcriptome-index /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.inx \
      --output-dir $BASEDIR/$sampleID \
      /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome $@
    if [ $? -eq 0 ]
    then
        echo "Finished Tophat mapping for $SampleID"
    else
        echo "Tophat mapping failed on $SampleID"
        exit 1
    fi    
    mkdir /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM
    mv $BASEDIR/$sampleID/accepted_hits.bam $BASEDIR/$sampleID/BAM/
    mv $BASEDIR/$sampleID/unmapped.bam $BASEDIR/$sampleID/BAM/
fi

if [ ! -f $BASEDIR/$sampleID/BAM/$sampleID.sort.bam ]
then
    echo "Sorting $sampleID"
    samtools sort $BASEDIR/$sampleID/BAM/accepted_hits.bam $BASEDIR/$sampleID/BAM/$sampleID.sort
    if [ $? -eq 0 ]
    then
        echo "Finished sorting for $SampleID"
    else
        echo "Could not sort BAM for $SampleID"
        exit 1
    fi
fi   

if [ ! -f $BASEDIR/$sampleID/BAM/$sampleID.sort.bam.bai ]
then
    echo "Indexing $sampleID"
    samtools index $BASEDIR/$sampleID/BAM/$sampleID.sort.bam    
    if [ $? -eq 0 ]
    then
        echo "Finished indexing for $SampleID"
    else
        echo "Could not index BAM for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$sampleID/BAM/$sampleID.chr.bam ]
then
    echo "Running RNAseqQC $sampleID"
    bash ~/LabNotes/SubmissionScripts/RNAseqQC.sh $BASEDIR/$sampleID/BAM/$sampleID.sort.bam   
    if [ $? -eq 0 ]
    then
        echo "Finished running RNAseqQC for $SampleID"
    else
        echo "Could not run RNAseqQC for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$sampleID/BAM/$sampleID.chr.counts.txt ]
then
    echo "Running htseq-count $sampleID"
    bash ~/LabNotes/SubmissionScripts/htseq-count.sh $BASEDIR/$sampleID/BAM/$sampleID.chr.bam   
    if [ $? -eq 0 ]
    then
        echo "Finished running htseq-count for $SampleID"
    else
        echo "Could not run htseq-count for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$sampleID/BAM/accepted_hits_fixup.bam ] || [ ! -f $BASEDIR/$sampleID/BAM/unmapped_fixup.bam ]
then
    echo "Fixing BAM formatting for $SampleID"
    bash ~/LabNotes/SubmissionScripts/tophat-recondition.sh $BASEDIR/$sampleID/BAM/
    if [ $? -eq 0 ]
    then
        echo "Finished BAM reformatting on $SampleID"
    else
        echo "tophat-recondition failed on $SampleID"
        exit 1
    fi    
fi

if [ ! -f $BASEDIR/$sampleID/BAM/accepted_hits_fixup_merge.bam ]
then
    echo "Merging mapped and unmapped BAM files for $SampleID"
    bash ~/LabNotes/SubmissionScripts/SamtoolsMerge.sh $BASEDIR/$sampleID/BAM/accepted_hits_fixup.bam $BASEDIR/$sampleID/BAM/unmapped_fixup.bam 
    if [ $? -eq 0 ]
    then
        echo "Finished BAM merging on $SampleID"
    else
        echo "Samtools merge failed on $SampleID"
        exit 1
    fi    
fi

if [ ! -f $BASEDIR/$sampleID/BAM/accepted_hits_fixup_merge_sort.bam ]
then
    echo "Sorting BAM for $SampleID"
    # Sort BAM files by query name
    bash ~/LabNotes/SubmissionScripts/SamtoolsSort.sh $BASEDIR/$sampleID/BAM/accepted_hits_fixup_merge.bam
    if [ $? -eq 0 ]
    then
        echo "Finished BAM sorting on $SampleID"
    else
        echo "Samtools sort failed on $SampleID"
        exit 1
    fi    
fi

if [ ! -f $BASEDIR/$sampleID/BAM/accepted_hits_fixup_merge_sort_RG.bam ]
then
    echo "Adding read groups for $SampleID"
    bash ~/LabNotes/SubmissionScripts/AddRG.sh $BASEDIR/$sampleID/BAM/accepted_hits_fixup_merge_sort.bam $SampleID
    if [ $? -eq 0 ]
    then
        echo "Finished adding read groups for $SampleID"
    else
        echo "could not add read groups for $SampleID"
        exit 1
    fi
fi


#bash ~/LabNotes/SubmissionScripts/dexseq-count.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.bam
#bash ~/LabNotes/SubmissionScripts/DivideBAM.sh $sampleID
#bash ~/LabNotes/SubmissionScripts/CallSNPs.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/Chromosomes/$sampleID.chr22.bam
#bash ~/LabNotes/SubmissionScripts/GTcheck.sh $sampleID
#index=`grep $sampleID ~/LabNotes/VCFindex.txt | cut -f 2`
#bash ~/LabNotes/SubmissionScripts/WASPnonRef.sh $sampleID $index
#bash ~/LabNotes/SubmissionScripts/RNAseqQCwasp.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.sorted.bam
#samtools sort /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.dedup.bam /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.dedup.sort
#bash ~/LabNotes/SubmissionScripts/RNAseqQCwasp.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.dedup.sort.bam
#bash ~/LabNotes/SubmissionScripts/clipOverlap.sh /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.dedup.sort.bam
#samtools index /c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM/$sampleID.chr.nonref.merged.dedup.sort.clip.bam
