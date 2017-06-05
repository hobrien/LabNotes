#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
for input in $@
#for input in `find /c8000xd3/rnaseq-heath/Mappings/ -name *.chr.nonref.merged.dedup.sort.clip.bam`
do
    
    echo "Started processing $input"
    python ~/src/WASP-0.2.1/mapping/rmdup_pe.py $input ${input%.*}_dedup.bam 
    echo "Finished processing $input"
done
exit $?
