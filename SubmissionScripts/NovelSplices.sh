#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH


java -jar ~/src/QoRTs.jar mergeNovelSplices  \
                --minCount 6 \
                --stranded \
                /c8000xd3/rnaseq-heath/Carolina/Counts/ \
                ~/LabNotes/Carolina_sampleInfo.txt \
                /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
                /c8000xd3/rnaseq-heath/Carolina/Counts/
                