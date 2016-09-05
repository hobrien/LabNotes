#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH
export PERL5LIB=/home/mpmho/perl5/lib/perl5:/usr/lib64/perl5
cd /c8000xd3/rnaseq-heath/DB
~/src/ncbi-blast-2.3.0+/bin/update_blastdb.pl nt
bash ~/LabNotes/SubmissionScripts/Untar.sh *.tar.gz
