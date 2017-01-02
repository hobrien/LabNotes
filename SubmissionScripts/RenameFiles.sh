#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
nameStem=$1
#mv $nameStem $nameStem-1
for filename in `find 12546-1 -name $nameStem.*`
do 
    echo ${filename#$nameStem}
done
