#Exeter
for id in 17068 17109 16929 #16840 
do
    find /c8000xd3/databank/foetal-rna/ -name $id*.fastq | xargs qsub ~/SubmissionScripts/Tophat2ex.sh
done

#Edinburgh
for id in  #18655 18294 18687 16438 16488 16286 18349 17049 18596 17087 17130 17701 17115 16972 17054 #17081 15768 15468 17812
do
    folder=`find /c8000xd3/databank/foetal-rna/ -name $id`
    find $folder -name *fastq.gz | xargs qsub ~/SubmissionScripts/Tophat2ed.sh
done
