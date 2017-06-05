#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=20G
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
for input in $@
#for input in `find /c8000xd3/rnaseq-heath/Mappings/ -name *.chr.nonref.merged.dedup.sort.clip.bam`
do
    
    echo "Started processing $input"
    python ~/src/WASP-0.2.1/mapping/find_intersecting_snps.py \
          --is_paired_end \
          --is_sorted \
          --output_dir ${input%/*}/find_intersecting_snps/ \
          --snp_tab /c8000xd3/rnaseq-heath/Genotypes/Imputation3/HDF5/snp_tab.h5 \
          --snp_index /c8000xd3/rnaseq-heath/Genotypes/Imputation3/HDF5/snp_index.h5 \
          --haplotype /c8000xd3/rnaseq-heath/Genotypes/Imputation3/HDF5/haplotypes.h5 \
          $input
    echo "Finished processing $input"
done
exit $?
