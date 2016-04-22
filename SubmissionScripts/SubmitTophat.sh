#Exeter
for id in #15641 16024 16115 16385 16428 16491 16548 16826 17048 17053 17071 17921 #15533 16929 16840 17068 17109  
do
  if [ ! -d /c8000xd3/rnaseq-heath/Mappings/$id ]
  then
    find /c8000xd3/databank/foetal-rna/ -name $id*.fastq | xargs qsub ~/SubmissionScripts/Tophat2ex.sh
  fi
done

#Compressed Exeter data
for id in #16810 
do
  if [ ! -d /c8000xd3/rnaseq-heath/Mappings/$id ]
  then
    find /c8000xd3/databank/foetal-rna/ -name $id*.fastq.gz | xargs qsub ~/SubmissionScripts/Tophat2ex.sh
  fi
done

#Edinburgh
for id in 17198 17229 17333 17475 #15655 16483 16640 17013 17025 #15468 15768 16286 16438 16488 16972 17049 18596 17701 17054 17081 17087 17115 17130 17812 18655 18294 18687 18349 
do
  if [ ! -d /c8000xd3/rnaseq-heath/Mappings/$id ]
  then
    folder=`find /c8000xd3/databank/foetal-rna/ -name $id`
    find $folder -name *fastq.gz | xargs qsub ~/SubmissionScripts/Tophat2ed.sh
  fi
done
