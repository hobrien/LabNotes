#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
trim_galore --illumina --paired --clip_R1 5 --clip_R2 5 --output_dir /c8000xd3/rnaseq-heath/Trimmed/ $@
