#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

cd /c8000xd3/rnaseq-heath/DB
ncbi-blast-2.3.0+/bin/update_blastdb.pl nt