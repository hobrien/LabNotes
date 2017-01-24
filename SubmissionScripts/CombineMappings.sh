#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=20G
#

BASEDIR=/c8000xd3/rnaseq-heath/Mappings/

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
SampleID=$1
#bash ~/LabNotes/SubmissionScripts/RenameFiles $SampleID
echo "Combining mappings for $SampleID"
if [ ! -d $BASEDIR/$SampleID/BAM ]
then
    echo "Making directory $BASEDIR/$SampleID/BAM"
    mkdir -p $BASEDIR/$SampleID/BAM
    if [ $? -eq 0 ]
    then
        echo "Finished making directory $BASEDIR/$SampleID/BAM"
    else
        echo "Could make directory $BASEDIR/$SampleID/BAM"
        exit 1
    fi
fi
if [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.chr.counts.txt ]
then
    echo "Combining counts for $SampleID"
    Rscript ~/LabNotes/R/CombineCounts.R $SampleID
    if [ $? -eq 0 ]
    then
        echo "Finished combining counts for $SampleID"
    else
        echo "Could combine counts for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.bam ]
then
    echo "Merging mappings for $SampleID"
    bash ~/LabNotes/SubmissionScripts/MergeSAM.sh `find $BASEDIR -name accepted_hits_filtered_sort_dedup_sort_RG.bam | grep $SampleID- | sort`

    if [ $? -eq 0 ]
    then
        echo "Finished merging mappings for $SampleID"
    else
        echo "Could merge mappings for $SampleID"
        exit 1
    fi
    mv $BASEDIR/$SampleID-1/BAM/accepted_hits_filtered_sort_dedup_sort_RG_merge.bam $BASEDIR/$SampleID/BAM/$SampleID.bam
fi

if [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.bam.bai ]
then
    echo "Indexing merged mappings for $SampleID"
    bash ~/LabNotes/SubmissionScripts/SamtoolsIndex.sh $BASEDIR/$SampleID/BAM/$SampleID.bam

    if [ $? -eq 0 ]
    then
        echo "Finished indexing merged mappings for $SampleID"
    else
        echo "Could index merged mappings for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/$SampleID.in.stats.txt ]
then
    echo "Indexing merged mappings for $SampleID"
    bash ~/LabNotes/SubmissionScripts/RNAseqQC.sh $BASEDIR/$SampleID/BAM/$SampleID.bam

    if [ $? -eq 0 ]
    then
        echo "Finished indexing merged mappings for $SampleID"
    else
        echo "Could index merged mappings for $SampleID"
        exit 1
    fi
fi

for chr in {1..22}
do
    if [ ! -f /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.filter_dup.vcf.gz ]
    then
        echo "Running ProcessVCF on $chr"
        bash ~/LabNotes/SubmissionScripts/ProcessVCF.sh $chr
        if [ $? -eq 0 ]
        then
            echo "Finished running ProcessVCF on $chr"
        else
            echo "Could not run ProcessVCF on $chr"
            exit 1
        fi
    fi
    if [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.$chr.ase.rtable ]
    then
        echo "Running ASEReadCounter on $chr for $SampleID"
        bash ~/LabNotes/SubmissionScripts/ASEReadCounter.sh $BASEDIR/$SampleID/BAM/$SampleID.bam $chr
        if [ $? -eq 0 ]
        then
            echo "Finished running ASEReadCounter on $chr for $SampleID"
        else
            echo "Could not run ASEReadCounter on $chr for $SampleID"
            exit 1
        fi
    fi
done
echo "Finished combining mappings for $SampleID"
exit $?

