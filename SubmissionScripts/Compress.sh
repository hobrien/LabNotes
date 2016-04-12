#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
filename=$(basename "$@")
gzip -c "$@" > /home/mpmho/Raw/${filename}.gz
