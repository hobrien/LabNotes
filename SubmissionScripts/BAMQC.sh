#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash

bamqc --outdir=$@ --gff /home/mpmho/Index/knownGene.gtf $@/accepted_hits.bam
