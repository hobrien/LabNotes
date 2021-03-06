#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

BASEDIR=/c8000xd3/rnaseq-heath/Genotypes/Imputation3
cd $BASEDIR/HDF5

echo "Running VCF2HDF5"
for chr in {1..22}
do
    if [ ! -f $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.vcf.gz ]
    then
        echo "Running ProcessVCF on $chr"
        bash ~/LabNotes/SubmissionScripts/ProcessVCF.sh $chr
        if [ $? -eq 0 ]
        then
            echo "Finished running ProcessVCF on $chr"
        else
            echo "Could not run ProcessVCF on $chr"
            exit 1
        fi
    fi
done


if [ ! -f $BASEDIR/HDF5/haplotypes.h5 ] | \
   [ ! -f $BASEDIR/HDF5/snp_index.h5 ] | \
   [ ! -f $BASEDIR/HDF5/snp_tab.h5 ]
then
    echo "Converting VCF files with duplicated sites to HDF5"
    ~/src/WASP-0.2.1/snp2h5/snp2h5 \
        --chrom /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/chromInfo.txt \
        --format vcf \
        --haplotype $BASEDIR/GRCh38/HDF5/haplotypes.h5 \
        --snp_index $BASEDIR/GRCh38/HDF5/snp_index.h5 \
        --snp_tab $BASEDIR/GRCh38/HDF5/snp_tab.h5 \
      $BASEDIR/GRCh38/chr*dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.vcf.gz
    if [ $? -eq 0 ]
    then
        echo "Finished converting VCF files with duplicated sites to HDF5"
    else
        echo "Could not convert VCF files with duplicated sites to HDF5"
        exit 1
    fi
fi
echo "Finished running VCF2HDF5"
exit $?
    
