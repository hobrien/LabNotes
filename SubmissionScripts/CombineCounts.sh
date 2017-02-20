#for SampleID in 11875 12107 12545 12546 12993 12994 13008 13142 15533 16483 17046 17054 17221 17264 18153 18540 18655 A138  A226  
LABNOTES='/Users/heo3/BTSync/FetalRNAseq/LabNotes'
for Individual in 1047 1102 1118 11305 11309 11349 11373 11423 11424 11429 11449 11451 11457 11489 11494 11496 11501 11511 11572 11654 11666 11683 11775 11794 11817 11834 11844 11854 11885 11892 11900 11947 12007 1290 1558 1835 1923 
do
    for SampleID in `grep "\b$Individual\b" $LABNOTES/sequences.txt | grep _1 | cut -f 1`
    do
        echo $SampleID
        grep $SampleID ${LABNOTES}/chr_bam_R1.txt | cut -f 1,2 | perl -pe 's/\s+/BAM\//' | perl -pe 's/bam/counts.txt/' | sort | uniq | xargs Rscript $LABNOTES/R/CombineCounts.R  ${LABNOTES}/Counts/raw/$SampleID
    done
done
