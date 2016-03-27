#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
#/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=/home/heath/Ref/hg19.fa O=/home/heath/Ref/hg19.dict

#/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar ReorderSam INPUT=/home/heath/Mappings/15533_300/accepted_hits.bam OUTPUT=/home/heath/Mappings/15533_300/accepted_hits_sorted.bam REFERENCE=/home/heath/Ref/hg19.fa

#/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar CollectRnaSeqMetrics REF_FLAT=/home/heath/Ref/refFlat.txt RIBOSOMAL_INTERVALS=/home/heath/Ref/hg19.rRNA.interval_list STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND INPUT=/home/heath/Mappings/15533_300/accepted_hits_sorted.bam OUTPUT=/home/heath/Mappings/15533_300/RnaSeqMetrics.txt ASSUME_SORTED=false

#bam_stat.py -i /home/heath/Mappings/15533_300/accepted_hits.bam > /home/heath/Mappings/15533_300/accepted_hits_stats.txt

# determine the strand of experiment ("1++,1--,2+-,2-+" = first strand, "1+-,1-+,2++,2--" = second strand)
#infer_experiment.py -r /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i /home/heath/Mappings/15533_300/accepted_hits.bam > /home/heath/Mappings/15533_300/accepted_hits_expt.txt

# plot distribution of insert sizes (size - total read length)
inner_distance.py -r /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i /home/heath/Mappings/15533_300_secondstrand/accepted_hits.bam -o /home/heath/Mappings/15533_300_secondstrand/15533_300_secondstrand -u 5000 -s 50 >/dev/null
/share/apps/R-3.2.2/bin/Rscript /home/heath/Mappings/15533_300_secondstrand/15533_300_secondstrand.inner_distance_plot.r

#junction_annotation.py -r /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i /home/heath/Mappings/15533_300_secondstrand/accepted_hits.bam -o /home/heath/Mappings/15533_300_secondstrand/15533_300_secondstrand
