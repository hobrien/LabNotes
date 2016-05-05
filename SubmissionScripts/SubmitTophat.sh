#Exeter
for id in 15533 #15533_2 15641 16024 16115 16385 16428 16491 16548 16826 17048 17053 17071 17921 # 16929 16840 17068 17109  
do
  if [ ! -d /c8000xd3/rnaseq-heath/Mappings/$id ]
  then
    find /c8000xd3/databank/foetal-rna/ -name $id*.fastq | xargs qsub ~/LabNotes/SubmissionScripts/Tophat2ex.sh
  fi
done

#Compressed Exeter data
for id in #16810 15533
do
  if [ ! -d /c8000xd3/rnaseq-heath/Mappings/$id ]
  then
    find /c8000xd3/databank/foetal-rna/ -name $id*.fastq.gz | xargs qsub ~/LabNotes/SubmissionScripts/Tophat2ex.sh
  fi
done

#Edinburgh
for id in #18432 18559 18694 18153 18241 18266 18282 18355 18134 17369 17671 17922 17923 18055 18121 A138 A226 15240 16649 17072 17160 17167 17175 #17629 17701 17812 17835 18249 18349 18372 18655 18666 18687 18983 19043 19052 #15655 16483 16640 17013 17025 17198 17229 17333 17475 17543 #15468 15768 16286 16438 16488 16972 17049 18596 17701 17054 17081 17087 17115 17130 17812 18655 18294 18687 18349 
do
  if [ ! -d /c8000xd3/rnaseq-heath/Mappings/$id ]
  then
    folder=`find /c8000xd3/databank/foetal-rna/ -name $id`
    find $folder -name *fastq.gz | xargs qsub ~/LabNotes/SubmissionScripts/Tophat2ed.sh
  fi
done
