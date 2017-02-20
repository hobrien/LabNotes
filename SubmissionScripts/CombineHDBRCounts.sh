LABNOTES='/Users/heo3/BTSync/FetalRNAseq/LabNotes'

for SampleID in 1047 1102 1118 11305 11309 11349 11373 11423 11424 11429 11449 11451 11457 11489 11494 11496 11501 11511 11572 11654 11666 11683 11775 11794 11817 11834 11844 11854 11885 11892 11900 11947 12007 1290 1558 1835 1923 
do
     grep "\b$SampleID\b" ${LABNOTES}/sequences.txt | cut -f 1 | grep _1 | perl -pe 's/1$/counts.txt/' | perl -pe 's/^/\/Users\/heo3\/BTSync\/FetalRNAseq\/Counts\/raw\//' | xargs Rscript $LABNOTES/R/CombineCounts.R /Users/heo3/BTSync/FetalRNAseq/Counts/raw/$SampleID.counts.txt
done
