#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

java -Xmx12g -jar ~/src/picard-tools-1.139/picard.jar ValidateSamFile I=$1 #MODE=SUMMARY
