#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
file=$@

if [[ $file == *"tar.gz" ]]
then
    tar -xzf $file
elif [[ $file == *"tar.bz2" ]]
then
    tar -xjf $file
elif [[ $file == *"tar" ]]
then
    tar -xf $file
else
    echo "extension not recognised"
fi   
 
