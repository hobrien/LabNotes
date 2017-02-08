for SampleID in 11875 12107 12545 12546 12993 12994 13008 13142 15533 16483 17046 17054 17221 17264 18153 18540 18655 A138  A226 
do
    grep $SampleID ~/LabNotes/chr_bam_R1.txt | cut -f 1,2 | perl -pe 's/\s+/BAM\//' | perl -pe 's/bam/counts.txt/' | sort | uniq | xargs Rscript ~/LabNotes/R/CombineCounts.R 
done
