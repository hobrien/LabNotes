#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export FONTCONFIG_PATH=/etc/fonts
Rscript ~/LabNotes/R/FvsM.r
