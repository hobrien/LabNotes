#Exeter
for id in 15641 16024 16115 16385 16428 16491 16548 16826 17048 17053 17071 17921 #15533 16929 16840 17068 17109  
do
  if [ ! -f /c8000xd3/rnaseq-heath/Mappings/$id/BAM/$id.sort.bam ]
  then
    find /c8000xd3/databank/foetal-rna/ -name $id*.fastq | xargs qsub ~/SubmissionScripts/Tophat2ex.sh
  fi
done

#Compressed Exeter data
for id in 16810 
do
  if [ ! -f /c8000xd3/rnaseq-heath/Mappings/$id/BAM/$id.sort.bam ]
  then
    find /c8000xd3/databank/foetal-rna/ -name $id*.fastq.gz | xargs qsub ~/SubmissionScripts/Tophat2ex.sh
  fi
done

#Edinburgh
for id in  #18655 18294 18687 16438 16488 16286 18349 17049 18596 17087 17130 17701 17115 16972 17054 #17081 15768 15468 17812
do
    folder=`find /c8000xd3/databank/foetal-rna/ -name $id`
    find $folder -name *fastq.gz | xargs qsub ~/SubmissionScripts/Tophat2ed.sh
done
