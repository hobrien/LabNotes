#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
for file in $@
do
    if [[ $file == *"tar.gz" ]]
    then
        tar -xzf $file
    elif [[ $file == *"tar.bz2" ]]
    then
        tar -xjf $file
    elif [[ $file == *"tar" ]]
    then
        tar -xf $file
    elif [[ $file == *"gz" ]]
    then
        gunzip $file
    elif [[ $file == *"zip" ]]
    then
        unzip $file
    else
        echo "extension not recognised"
    fi   
done 
