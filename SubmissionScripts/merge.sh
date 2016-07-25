#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

find /c8000xd3/rnaseq-heath/Mappings/ -maxdepth 3 -name *chr.nonref.merged.dedup.sort.bam | xargs samtools merge -r /c8000xd3/rnaseq-heath/Mappings/All.bam $@
