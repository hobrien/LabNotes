# command to extract from all samples: find /c8000xd3/rnaseq-heath/Mappings/ -maxdepth 2 -name junctions.bed | grep -v SRR | grep -v 50bp | grep -v _ | xargs bash SubmissionScripts/ExtractHTTjuncitons.bed


for dataset in $@
do
    folder_path=${dataset%/*}
    sampleID=${folder_path##*/}

    echo $sampleID
done
