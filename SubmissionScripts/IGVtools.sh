#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
filename1=${1##*/}
sampleID=${filename1%%.*}


/share/apps/IGV_2.3.72/IGVTools/igvtools sort /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/$sampleID.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/$sampleID.sorted.gtf

/share/apps/IGV_2.3.72/IGVTools/igvtools index /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/$sampleID.sorted.gtf
