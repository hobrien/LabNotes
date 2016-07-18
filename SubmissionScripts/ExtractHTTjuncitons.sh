# command to extract from all samples: find /c8000xd3/rnaseq-heath/Mappings/ -maxdepth 2 -name junctions.bed | grep -v SRR | grep -v 50bp | grep -v _ | xargs bash SubmissionScripts/ExtractHTTjuncitons.bed


for dataset in $@
do
    folder_path=${dataset%/*}
    sampleID=${folder_path##*/}

    grep chr4 /c8000xd3/rnaseq-heath/Mappings/$sampleID/junctions.bed \
      | awk '{ if ($2 > 3074680 && $3 < 3243960) print $0}' \
      > ~/HTT/$sampleID.bed
done