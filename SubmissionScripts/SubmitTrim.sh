#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
if [ ! -f  ~/Trimmed/${1##*/}_val_1.fq.gz -a ! -f ~/Trimmed/${2##*/}_val_2.fq.gz ]
then
    qsub ~/SubmissionScripts/Trim.sh $@
fi
