#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
gtf2bed.pl -x $@ > ${@%.*}.temp.bed 
perl -ne 'print if /chr\d\t/' ${@%.*}.temp.bed > ${@%.*}.bed
perl -ne 'print if /chr\d\d\t/' ${@%.*}.temp.bed >> ${@%.*}.bed
perl -ne 'print if /chr[MXY]\t/' ${@%.*}.temp.bed >> ${@%.*}.bed
perl -ne 'print if /chr\d_/' ${@%.*}.temp.bed >> ${@%.*}.bed
perl -ne 'print if /chr\d\d_/' ${@%.*}.temp.bed >> ${@%.*}.bed
perl -ne 'print if /chr[MXY]_/' ${@%.*}.temp.bed >> ${@%.*}.bed
perl -ne 'print if /chrUn_/' ${@%.*}.temp.bed >> ${@%.*}.bed
perl -ne 'print if not /chr/' ${@%.*}.temp.bed >> ${@%.*}.bed
rm ${@%.*}.temp.bed