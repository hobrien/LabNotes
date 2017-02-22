#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
for file in $@
do
    echo "Downloading $file"
    curl -O -L "$file"
    if [ $? -eq 0 ]
    then
        echo "Finished downloading $file"
    else
        echo "Could not download $file"
        exit 1
    fi    
done 
