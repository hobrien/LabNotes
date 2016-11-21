sampleID=$1
for num in {1..22}
do
    qsub ~/LabNotes/SubmissionScripts/WASPnonRefchr.sh $sampleID $num
done
