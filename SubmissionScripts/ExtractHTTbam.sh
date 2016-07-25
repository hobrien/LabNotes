# command to extract from all samples: find /c8000xd3/rnaseq-heath/Mappings/ -maxdepth 3 -name *.chr.bam | grep -v SRR | grep -v 50bp | grep -v _ | xargs bash SubmissionScripts/ExtractHTTjuncitons.bed


for dataset in $@
do
    folder_path=${dataset%/*}
    folder_path=${folder_path%/*}
    sampleID=${folder_path##*/}

    samtools view -b $dataset chr4:3074680-3243960 >~/HTT/$sampleID.htt.bam
done
