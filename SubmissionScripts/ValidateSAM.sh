#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

java -Xmx2g -jar ~/src/picard-tools-1.139/picard.jar ValidateSamFile I=/c8000xd3/rnaseq-heath/Mappings/15240/BAM/15240.chr.nonref.merged.dedup.bam MODE=SUMMARY
