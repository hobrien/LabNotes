#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

java -jar ~/src/picard-tools-1.139/picard.jar FixMateInformation \
  I=$1 \
  O=fixed_mate_{$1}
