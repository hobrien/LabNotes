#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

python ~/src/tophat-recondition/tophat-recondition.py -q $1
exit $?
