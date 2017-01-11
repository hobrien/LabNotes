#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

BASEDIR=/c8000xd3/rnaseq-heath/Mappings

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
path=${1%/*}
path=${path%/*}
SampleID=${path##*/}
filename=${1##*/}
echo "Starting WASP Remapping on $SampleID"

mkdir $BASEDIR/$SampleID/BAM/find_intersecting_snps

if [ ! -f $BASEDIR/$SampleID/BAM/find_intersecting_snps/${filename%.*}.remap.fq1.gz ] || [ ! -f $BASEDIR/$SampleID/BAM/find_intersecting_snps/${filename%.*}.remap.fq2.gz ]
then
    echo "finding intersecting SNPs for $SampleID"
    python ~/src/WASP-0.2.1/mapping/find_intersecting_snps.py \
          --is_paired_end \
          --is_sorted \
          --output_dir $BASEDIR/$SampleID/BAM/find_intersecting_snps/ \
          --snp_tab /c8000xd3/rnaseq-heath/Genotypes/Imputation3/HDF5/snp_tab.h5 \
          --snp_index /c8000xd3/rnaseq-heath/Genotypes/Imputation3/HDF5/snp_index.h5 \
          --haplotype /c8000xd3/rnaseq-heath/Genotypes/Imputation3/HDF5/haplotypes.h5 \
          $1
    if [ $? -eq 0 ]
    then
        echo "Finished finding intersecting SNPs for $SampleID"
    else
        echo "Could not find intersecting SNPs $SampleID"
        exit 1
    fi  
fi  

if [ ! -f $BASEDIR/$SampleID/BAM/remap_intersecting_snps/accepted_hits.bam ]
then
    echo "remapping reads with intersecting SNPs for $SampleID"
    mkdir $BASEDIR/$SampleID/BAM/remap_intersecting_snps
    tophat --keep-fasta-order --library-type fr-secondstrand --mate-inner-dist 500 --mate-std-dev 50 --num-threads 8 \
          --transcriptome-index /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.inx \
          --output-dir $BASEDIR/$SampleID/BAM/remap_intersecting_snps \
          /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome \
          $BASEDIR/$SampleID/BAM/find_intersecting_snps/${filename%.*}.remap.fq1.gz \
          $BASEDIR/$SampleID/BAM/find_intersecting_snps/${filename%.*}.remap.fq2.gz

    if [ $? -eq 0 ]
    then
        echo "Finished remapping reads for $SampleID"
    else
        echo "Could not remap reads for $SampleID"
        exit 1
    fi
fi
     
if [ ! -f $BASEDIR/$SampleID/BAM/remap_intersecting_snps/accepted_hits_sort.bam ]
then
    echo "Sorting remapped BAM $SampleID"
    bash ~/LabNotes/SubmissionScripts/SamtoolsSort.sh $BASEDIR/$SampleID/BAM/remap_intersecting_snps/accepted_hits.bam
    if [ $? -eq 0 ]
    then
        echo "Finished sorting remapped BAM for $SampleID"
    else
        echo "Could not sort remapped BAM for $SampleID"
        exit 1
    fi
fi   

if [ ! -f $BASEDIR/$SampleID/BAM/remap_intersecting_snps/accepted_hits_sort.bam.bai ]
then
    echo "Indexing remapped BAM for $SampleID"
    samtools index $BASEDIR/$SampleID/BAM/remap_intersecting_snps/accepted_hits_sort.bam    
    if [ $? -eq 0 ]
    then
        echo "Finished indexing remapped BAM for $SampleID"
    else
        echo "Could not index remapped BAM for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/BAM/remap_intersecting_snps/accepted_hits_sort_keep.bam ]
then
    echo "Filtering remapped reads for $SampleID"
    python ~/src/WASP-0.2.1/mapping/filter_remapped_reads.py \
          $BASEDIR/$SampleID/BAM/find_intersecting_snps/${filename%.*}.to.remap.bam \
          $BASEDIR/$SampleID/BAM/remap_intersecting_snps/accepted_hits_sort.bam \
          $BASEDIR/$SampleID/BAM/remap_intersecting_snps/accepted_hits_sort_keep.bam
    if [ $? -eq 0 ]
    then
        echo "Finished filtering remapped reads for $SampleID"
    else
        echo "Could not filter remapped reads for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/BAM/${filename%.*}_filtered.bam ]
then
    echo "Merging filtered reads for $SampleID"
    samtools merge $BASEDIR/$SampleID/BAM/${filename%.*}_filtered.bam \
              $BASEDIR/$SampleID/BAM/remap_intersecting_snps/accepted_hits_sort_keep.bam  \
              $BASEDIR/$SampleID/BAM/find_intersecting_snps/${filename%.*}.keep.bam
    if [ $? -eq 0 ]
    then
        echo "Finished merging filtered reads for $SampleID"
    else
        echo "Could not merge filtered reads for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/BAM/${filename%.*}_filtered_sort.bam ]
then
    echo "Sorting filtered BAM $SampleID"
    bash ~/LabNotes/SubmissionScripts/SamtoolsSort.sh $BASEDIR/$SampleID/BAM/${filename%.*}_filtered.bam
    if [ $? -eq 0 ]
    then
        echo "Finished sorting filtered BAM for $SampleID"
    else
        echo "Could not sort filtered BAM for $SampleID"
        exit 1
    fi
fi   

if [ ! -f $BASEDIR/$SampleID/BAM/${filename%.*}_filtered_sort.bam.bai ]
then
    echo "Indexing filtered BAM for $SampleID"
    samtools index $BASEDIR/$SampleID/BAM/${filename%.*}_filtered_sort.bam  
    if [ $? -eq 0 ]
    then
        echo "Finished indexing filtered BAM for $SampleID"
    else
        echo "Could not index filtered BAM for $SampleID"
        exit 1
    fi
fi

if [ ! -f $BASEDIR/$SampleID/BAM/${filename%.*}_filtered_sort_dedup.bam ]
then
    echo "Deduplicating filtered BAM for $SampleID"
    python ~/src/WASP-0.2.1/mapping/rmdup_pe.py \
        $BASEDIR/$SampleID/BAM/${filename%.*}_filtered_sort.bam  \
        $BASEDIR/$SampleID/BAM/${filename%.*}_filtered_sort_dedup.bam
    if [ $? -eq 0 ]
    then
        echo "Finished deduplicating filtered BAM for $SampleID"
    else
        echo "Could not deduplicate filtered BAM for $SampleID"
        exit 1
    fi
fi

