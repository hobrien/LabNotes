#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:$PATH
folder_path=${@%/*}
folder_path=${folder_path%/*}
folder=${folder_path##*/}

# [BAMQC](https://github.com/s-andrews/BamQC)
#bamqc --outdir=$folder_path --gff /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf $@

# [PicardTools](http://broadinstitute.github.io/picard/)
#/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=/home/heath/Ref/hg19.fa O=/home/heath/Ref/hg19.dict
#/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar ReorderSam INPUT=/home/heath/Mappings/15533_300/accepted_hits.bam OUTPUT=/home/heath/Mappings/15533_300/accepted_hits_sorted.bam REFERENCE=/home/heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa
#/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar CollectRnaSeqMetrics REF_FLAT=/home/heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/refFlat.txt.gz STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND INPUT=/home/heath/Mappings/15533_300_secondstrand/accepted_hits.bam OUTPUT=/home/heath/Mappings/15533_300_secondstrand/RnaSeqMetrics.txt ASSUME_SORTED=false

# This isn't needed because it's run separately for reads mapping and not mapping to rDNA
bam_stat.py -i $@ > $folder_path/$folder.stats.txt

# determine the strand of experiment ("1++,1--,2+-,2-+" = first strand, "1+-,1-+,2++,2--" = second strand)
infer_experiment.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $@ > $folder_path/$folder.expt.txt

# plot distribution of insert sizes (size - total read length)
inner_distance.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $@ -o $folder_path/$folder -u 1000 -s 10 >/dev/null

# the necessary output from this is going to the log file, not to $folder.junction.txt
junction_annotation.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $@ -o $folder_path/$folder  >  $folder_path/$folder.junction.txt

junction_saturation.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $@ -o $folder_path/$folder

read_distribution.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $@ > $folder_path/$folder.dist.txt

read_duplication.py -i $@ -o $folder_path/$folder

geneBody_coverage.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $@ -o $folder_path/$folder

split_bam.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/AbundantSequences/humRibosomal.bed -i $@ -o $folder_path/BAM/$folder
bam_stat.py -i $folder_path/BAM/$folder.in.bam > $folder_path/$folder.in.stats.txt
bam_stat.py -i $folder_path/BAM/$folder.ex.bam > $folder_path/$folder.ex.stats.txt

deletion_profile.py -l 120 -i $@ -o $folder_path/$folder
insertion_profile.py -s PE -i $@ -o $folder_path/$folder
#clipping_profile.py -s PE -i $@ -o $folder_path/$folder
mismatch_profile.py -l 120 -i $@ -o $folder_path/$folder

echo "finished QC for $@"
