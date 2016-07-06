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
        folder_path=${folder_path}/Wasp_nonref
        folder=${folder%%.*}.nrwasp        
    else    
        folder_path=${folder_path}/Wasp
        folder=${folder%%.*}.wasp
    fi
    echo "Starting QC for $dataset"

    mkdir $folder_path


    geneBody_coverage.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $dataset -o $folder_path/$folder



    echo "finished QC for $dataset"
done
