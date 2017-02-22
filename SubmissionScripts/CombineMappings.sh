#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=20G
#

BASEDIR=/c8000xd3/rnaseq-heath/ASEmappings/

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
SampleID=$1
#bash ~/LabNotes/SubmissionScripts/RenameFiles $SampleID
echo "Combining mappings for $SampleID"
if [ ! -d $BASEDIR/$SampleID ]
then
    echo "Making directory $BASEDIR/$SampleID"
    mkdir $BASEDIR/$SampleID
    if [ $? -eq 0 ]
    then
        echo "Finished making directory $BASEDIR/$SampleID"
    else
        echo "Could make directory $BASEDIR/$SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/$SampleID.bam ]
then
    echo "Merging mappings for $SampleID"
    bash ~/LabNotes/SubmissionScripts/MergeSAM.sh `find $BASEDIR -name ${SampleID}-*_filtered_dedup_sort_RG.bam | grep $SampleID- | sort`

    if [ $? -eq 0 ]
    then
        echo "Finished merging mappings for $SampleID"
    else
        echo "Could merge mappings for $SampleID"
        exit 1
    fi
    mv $BASEDIR/$SampleID-1/${SampleID}-1_filtered_dedup_sort_RG_merge.bam $BASEDIR/$SampleID/$SampleID.bam
fi

if [ ! -f $BASEDIR/$SampleID/$SampleID.bam.bai ]
then
    echo "Indexing merged mappings for $SampleID"
    bash ~/LabNotes/SubmissionScripts/SamtoolsIndex.sh $BASEDIR/$SampleID/$SampleID.bam

    if [ $? -eq 0 ]
    then
        echo "Finished indexing merged mappings for $SampleID"
    else
        echo "Could index merged mappings for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/$SampleID_stats.txt ]
then
    echo "Indexing merged mappings for $SampleID"
    bam_stat.py -i $BASEDIR/$SampleID/${SampleID}.bam > $BASEDIR/$SampleID/${SampleID}_stats.txt

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
    if [ ! -f $BASEDIR/$SampleID/$SampleID.$chr.ase.rtable ]
    then
        echo "Running ASEReadCounter on $chr for $SampleID"
        bash ~/LabNotes/SubmissionScripts/ASEReadCounter.sh $BASEDIR/$SampleID/$SampleID.bam $chr
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

