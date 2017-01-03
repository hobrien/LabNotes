#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
nameStem=$1
mv /c8000xd3/rnaseq-heath/Mappings/$nameStem /c8000xd3/rnaseq-heath/Mappings/$nameStem-1
for filename in `find /c8000xd3/rnaseq-heath/Mappings/$nameStem-1 -name $nameStem.*`
do 
    pathname=${filename%/*}
    mv $filename $pathname/$nameStem-1.${filename#*.}
done
