#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:$PATH

for dataset in $@
do
    folder_path=${dataset%/*}
    folder=${dataset##*/}
    if [[ $dataset == *"nonref"* ]]
    then
        if [[ $dataset == *"dedup"* ]]
        then
            folder_path=${folder_path}/Wasp_dedup
            folder=${folder%%.*}.dedup
        else
            folder_path=${folder_path}/Wasp_nonref
            folder=${folder%%.*}.nrwasp
        fi                             
    else    
        folder_path=${folder_path}/Wasp
        folder=${folder%%.*}.wasp
    fi
    if [[ ! -e $dataset.bai ]]
    then
        echo "indexing $dataset"
        samtools index $dataset
        if [ $? -ne 0 ]
        then
            echo "sorting and indexing $dataset"
            nosort=${dataset/.sort}
            mv $dataset $nosort
            base=${nosort%%.*}
            samtools sort $nosort $base.sort
            dataset=$base.sort.bam
            samtools index $dataset
        fi    
    fi    
    echo "Starting QC for $dataset"

    mkdir $folder_path

    bam_stat.py -i $dataset > $folder_path/$folder.stats.txt

    # determine the strand of experiment ("1++,1--,2+-,2-+" = first strand, "1+-,1-+,2++,2--" = second strand)
    infer_experiment.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $dataset > $folder_path/$folder.expt.txt

    # plot distribution of insert sizes (size - total read length)
    inner_distance.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $dataset -o $folder_path/$folder -u 1000 -s 10 >/dev/null

    # the necessary output from this is going to the log file, not to $folder.junction.txt
    junction_annotation.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $dataset -o $folder_path/$folder  >  $folder_path/$folder.junction.txt

    junction_saturation.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $dataset -o $folder_path/$folder

    read_distribution.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $dataset > $folder_path/$folder.dist.txt

    geneBody_coverage.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $dataset -o $folder_path/$folder



    echo "finished QC for $dataset"
done
