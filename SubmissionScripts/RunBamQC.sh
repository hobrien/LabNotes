#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

cd ~/LabNotes/R
Rscript -e "library(knitr); rmarkdown::render('BamQCwasp.rmd')"
