#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH


# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
path=${1%/*}
path=${path%/*}
sampleID=${path##*/}

echo "Starting WASP Remapping on $sampleID"


python ~/src/WASP/mapping/find_intersecting_snps.py \
          --is_paired_end \
          --is_sorted \
          --output_dir find_intersecting_snps \
          --snp_tab /c8000xd3/rnaseq-heath/Genotypes/Imputation3/HDF5/snp_tab.h5 \
          --snp_index /c8000xd3/rnaseq-heath/Genotypes/Imputation3/HDF5/snp_index.h5 \
          --haplotype /c8000xd3/rnaseq-heath/Genotypes/Imputation3/HDF5/haplotype.h5 \
          --samples my_samples.txt \
          $1
