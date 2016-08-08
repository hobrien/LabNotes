#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

for sampleID in 15240 15468 15533 15533_2 15641 15655 15768 16024 16115 16286 16385 16428 16438 16483 16488 16491 16548 16640 16649 16810 16826 16840 16929 16972 17013 17025 17048 17049 17053 17054 17068 17071 17072 17081 17087 17109 17115 17130 17160 17167 17175 17198 17229 17333 17369 17475 17543 17629 17671 17701 17812 17835 17921-l1 17922 17923 18055 18121 18134 18153 18241 18249 18266 18282 18294 18349 18355 18372 18432 18559 18596 18655 18666 18687 18694 18983 19043 19052 A138 A226 
do
    folder_path=/c8000xd3/rnaseq-heath/Mappings/$sampleID/BAM
    chr=$1
    mkdir $folder_path/Chromosomes
    #echo $folder_path/Chromosomes/$sampleID.chr$chr.bam
    samtools view -bh $folder_path/$sampleID.sort.bam chr$chr > $folder_path/Chromosomes/$sampleID.chr$chr.bam
    samtools index $folder_path/Chromosomes/$sampleID.chr$chr.bam
done
