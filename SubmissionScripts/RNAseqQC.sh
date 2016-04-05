#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:$PATH
#/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=/home/heath/Ref/hg19.fa O=/home/heath/Ref/hg19.dict
folder_path=${@%/*}
folder=${folder_path##*/}

#bamqc --outdir=$folder_path --gff /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf $@

#/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar ReorderSam INPUT=/home/heath/Mappings/15533_300/accepted_hits.bam OUTPUT=/home/heath/Mappings/15533_300/accepted_hits_sorted.bam REFERENCE=/home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa

#/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar CollectRnaSeqMetrics REF_FLAT=/home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/refFlat.txt.gz STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND INPUT=/home/heath/Mappings/15533_300_secondstrand/accepted_hits.bam OUTPUT=/home/heath/Mappings/15533_300_secondstrand/RnaSeqMetrics.txt ASSUME_SORTED=false
#bam_stat.py -i $@ > $folder_path/$folder.stats.txt

# determine the strand of experiment ("1++,1--,2+-,2-+" = first strand, "1+-,1-+,2++,2--" = second strand)
#infer_experiment.py -r /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $@ > $folder_path/$folder.expt.txt

# plot distribution of insert sizes (size - total read length)
#inner_distance.py -r /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $@ -o $folder_path/$folder -u 1000 -s 10 >/dev/null

#junction_annotation.py -r /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $@ -o $folder_path/$folder

#junction_saturation.py -r /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $@ -o $folder_path/$folder

#read_distribution.py -r /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i /home/heath/Mappings/15533_300_secondstrand/accepted_hits.bam > -i $@ > $folder_path/$folder.dist.txt

#read_duplication.py -i $@ -o $folder_path/$folder

#geneBody_coverage.py -r /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $@ -o $folder_path/$folder

split_bam.py -r /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/AbundantSequences/humRibosomal.bed -i $@ -o $folder_path/$folder
bam_stat.py -i $folder_path/$folder.in.bam > $folder_path/$folder.in.stats.txt
bam_stat.py -i $folder_path/$folder.ex.bam > $folder_path/$folder.ex.stats.txt

