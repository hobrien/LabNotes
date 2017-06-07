#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution

# unique + non-unique + unmapped + non-primary = total records
# unique + non-unique + unmapped = total reads
# Read-1 + Read-2 = unique
# Reads map to '+' + Reads map to '-' = unique
for input in $@
#for input in `find /c8000xd3/rnaseq-heath/Mappings/ -name *.chr.nonref.merged.dedup.sort.clip.bam`
do
    
    echo "Started processing $input"
    bam_stat.py -i $input > ${input%.*}_stats.txt
    echo "Finished processing $input"
done
exit $?
