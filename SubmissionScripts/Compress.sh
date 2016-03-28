#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
filename=$(basename "$@")
gzip -c "$@" > /home/heath/Raw/${filename}.gz
