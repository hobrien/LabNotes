#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
curl -O -L "$@" 
