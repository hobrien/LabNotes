#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 8
#

# This will take an individual pair of read files and take them through all steps of
# mapping, remapping, deduplicating and read group addition.
# Once this has been run on all datasets for a sample, they can be combined and used for
# allele-specific read counting using CombineMappings.sh

#export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
SampleID=$1 
MAPPER=HISAT2
BASEDIR=/c8000xd3/rnaseq-heath/ASEmappings
ANNOTATION=/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode
REF=/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta

echo "Starting mapping for $BASEDIR/$SampleID"
if [ ! -d $BASEDIR/$SampleID ]
then
    echo "Creating directory $BASEDIR/$SampleID"
    mkdir $BASEDIR/$SampleID
    if [ $? -eq 0 ]
    then
        echo "Finished creating directory for $BASEDIR/$SampleID"
    else
        echo "Could not create directory for $BASEDIR/$SampleID"
        exit 1
    fi    
fi

if [ ! -f $BASEDIR/$SampleID/$SampleID.bam ]
then
    echo "$MAPPER mapping for $SampleID"
    sequences=$(for name in `grep -P "\s$SampleID(\s|$)"  ~/LabNotes/sequences.txt | cut -f 1`; do find /c8000xd2/foetalRNAseq/ /c8000xd3/databank/foetal-rna/ -name $name*f*q.gz; done)
    if [ ! $sequences ]
    then
        echo "trying without gzip"
        sequences=$(for name in `grep -P "\s$SampleID(\s|$)"  ~/LabNotes/sequences.txt | cut -f 1`; do find /c8000xd2/foetalRNAseq/ /c8000xd3/databank/foetal-rna/ -name $name*f*q; done)
    fi
    if [ ! $sequences ]
    then
        echo "could not find sequence files for $SampleID"
        exit 1
    fi
    echo "Read files: $sequences"
    sequences="$sequences $BASEDIR/$SampleID/$SampleID.bam" # Add output file to arguments
    qsub -N h${SampleID}_map ~/LabNotes/SubmissionScripts/HISAT2.sh $sequences
    if [ $? -eq 0 ]
    then
        echo "Finished $MAPPER mapping for $SampleID"
    else
        echo "$MAPPER mapping failed on $SampleID"
        rm $BASEDIR/$SampleID/$SampleID.bam
        exit 1
    fi    
fi

if [ ! -f $BASEDIR/$SampleID/${SampleID}_sort.bam ]
then
    echo "Sorting $SampleID"
    samtools sort -o $BASEDIR/$SampleID/${SampleID}_sort.bam $BASEDIR/$SampleID/$SampleID.bam
    if [ $? -eq 0 ]
    then
        echo "Finished sorting for $SampleID"
    else
        echo "Could not sort BAM for $SampleID"
        exit 1
    fi
fi   

if [ ! -f $BASEDIR/$SampleID/${SampleID}_sort.bam.bai ]
then
    echo "Indexing $SampleID"
    samtools index $BASEDIR/$SampleID/${SampleID}_sort.bam    
    if [ $? -eq 0 ]
    then
        echo "Finished indexing for $SampleID"
    else
        echo "Could not index BAM for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/${SampleID}_stats.txt ]
then
    echo "Running RNAseqQC $SampleID"
    bam_stat.py -i $BASEDIR/$SampleID/${SampleID}_sort.bam > $BASEDIR/$SampleID/${SampleID}_stats.txt
    if [ $? -eq 0 ]
    then
        echo "Finished running RNAseqQC for $SampleID"
    else
        echo "Could not run RNAseqQC for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_sort.remap.fq1.gz ] || [ ! -f $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_sort.remap.fq2.gz ]
then
    echo "finding intersecting SNPs for $SampleID"
    mkdir $BASEDIR/$SampleID/find_intersecting_snps
    python ~/src/WASP-0.2.1/mapping/find_intersecting_snps.py \
          --is_paired_end \
          --is_sorted \
          --output_dir $BASEDIR/$SampleID/find_intersecting_snps/ \
          --snp_tab /c8000xd3/rnaseq-heath/Genotypes/Imputation3/HDF5/snp_tab.h5 \
          --snp_index /c8000xd3/rnaseq-heath/Genotypes/Imputation3/HDF5/snp_index.h5 \
          --haplotype /c8000xd3/rnaseq-heath/Genotypes/Imputation3/HDF5/haplotypes.h5 \
          $BASEDIR/$SampleID/${SampleID}_sort.bam
    if [ $? -eq 0 ]
    then
        echo "Finished finding intersecting SNPs for $SampleID"
    else
        echo "Could not find intersecting SNPs $SampleID"
        exit 1
    fi  
fi  

if [ ! -f $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_remap.bam ]
then
    echo "remapping reads with intersecting SNPs for $SampleID"
    bash ~/LabNotes/SubmissionScripts/HISAT2.sh \
      $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_sort.remap.fq1.gz \
      $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_sort.remap.fq2.gz \
      | samtools view -S -bo $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_remap.bam -

    if [ $? -eq 0 ]
    then
        echo "Finished remapping reads for $SampleID"
    else
        echo "Could not remap reads for $SampleID"
        exit 1
    fi
fi
     
if [ ! -f $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_remap_sort.bam ]
then
    echo "Sorting remapped BAM $SampleID"
    samtools sort -o $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_remap_sort.bam \
      $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_remap.bam
      
    if [ $? -eq 0 ]
    then
        echo "Finished sorting remapped BAM for $SampleID"
    else
        echo "Could not sort remapped BAM for $SampleID"
        exit 1
    fi
fi   

if [ ! -f $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_remap_sort.bam.bai ]
then
    echo "Indexing remapped BAM for $SampleID"
    samtools index $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_remap_sort.bam    
    if [ $? -eq 0 ]
    then
        echo "Finished indexing remapped BAM for $SampleID"
    else
        echo "Could not index remapped BAM for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_remap_sort_keep.bam ]
then
    echo "Filtering remapped reads for $SampleID"
    python ~/src/WASP-0.2.1/mapping/filter_remapped_reads.py \
          $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_sort.to.remap.bam \
          $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_remap_sort.bam \
          $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_remap_sort_keep.bam 
    if [ $? -eq 0 ]
    then
        echo "Finished filtering remapped reads for $SampleID"
    else
        echo "Could not filter remapped reads for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/${SampleID}_sort_filtered.bam ]
then
    echo "Merging filtered reads for $SampleID"
    samtools merge $BASEDIR/$SampleID/${SampleID}_sort_filtered.bam \
              $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_remap_sort_keep.bam  \
              $BASEDIR/$SampleID/find_intersecting_snps/${SampleID}_sort.keep.bam
    if [ $? -eq 0 ]
    then
        echo "Finished merging filtered reads for $SampleID"
    else
        echo "Could not merge filtered reads for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/${SampleID}_filtered_sort.bam ]
then
    echo "Sorting filtered BAM $SampleID"
    samtools sort -o $BASEDIR/$SampleID/${SampleID}_filtered_sort.bam \
      $BASEDIR/$SampleID/${SampleID}_sort_filtered.bam 
    if [ $? -eq 0 ]
    then
        echo "Finished sorting filtered BAM for $SampleID"
    else
        echo "Could not sort filtered BAM for $SampleID"
        exit 1
    fi
fi   

if [ ! -f $BASEDIR/$SampleID/${SampleID}_filtered_sort.bam.bai ]
then
    echo "Indexing filtered BAM for $SampleID"
    samtools index $BASEDIR/$SampleID/${SampleID}_filtered_sort.bam 
    if [ $? -eq 0 ]
    then
        echo "Finished indexing filtered BAM for $SampleID"
    else
        echo "Could not index filtered BAM for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/${SampleID}_filtered_sort_dedup.bam ]
then
    echo "Deduplicating filtered BAM for $SampleID"
    python ~/src/WASP-0.2.1/mapping/rmdup_pe.py \
        $BASEDIR/$SampleID/${SampleID}_filtered_sort.bam  \
        $BASEDIR/$SampleID/${SampleID}_filtered_sort_dedup.bam
    if [ $? -eq 0 ]
    then
        echo "Finished deduplicating filtered BAM for $SampleID"
    else
        echo "Could not deduplicate filtered BAM for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/${SampleID}_filtered_dedup_sort.bam ]
then
    echo "Sorting filtered, deduplicated BAM for $SampleID"
    # Sort BAM files by query name
    samtools sort -o $BASEDIR/$SampleID/${SampleID}_filtered_dedup_sort.bam \
      $BASEDIR/$SampleID/${SampleID}_filtered_sort_dedup.bam 
      
    if [ $? -eq 0 ]
    then
        echo "Finished sortting filtered, deduplicated BAM for $SampleID"
    else
        echo "Could not sort filtered, deduplicated BAM for $SampleID"
        exit 1
    fi    
fi

if [ ! -f $BASEDIR/$SampleID/${SampleID}_filtered_dedup_sort.bam.bai ]
then
    echo "Indexing filtered, deduplicated BAM for $SampleID"
    samtools index $BASEDIR/$SampleID/${SampleID}_filtered_dedup_sort.bam 
    if [ $? -eq 0 ]
    then
        echo "Finished indexing filtered, deduplicated BAM for $SampleID"
    else
        echo "Could not index filtered, deduplicated BAM for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/${SampleID}_filtered_dedup_sort_RG.bam ]
then
    echo "Adding read groups for $SampleID"
    bash ~/LabNotes/SubmissionScripts/AddRG.sh $BASEDIR/$SampleID/${SampleID}_filtered_dedup_sort.bam $SampleID
    if [ $? -eq 0 ]
    then
        echo "Finished adding read groups for $SampleID"
    else
        echo "could not add read groups for $SampleID"
        exit 1
    fi
fi
echo "Finished mapping for $SampleID"
exit $?

