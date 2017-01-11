#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

BASEDIR=/c8000xd3/rnaseq-heath/Genotypes/Imputation3
ls /c8000xd3/rnaseq-heath/Mappings | cut -d- -f 1 | sort | uniq > ~/LabNotes/mappings.txt
for chr in {1..22}
do 
    echo "Modifying header of chromosome $chr"
    head -13 ~/LabNotes/header.txt > header_temp.txt
    echo "##contig=<ID=$chr>" >> header_temp.txt 
    tail -1 ~/LabNotes/header.txt >> header_temp.txt
    bcftools reheader -h header_temp.txt \
        -o $BASEDIR/hg19/chr$chr.dose.rename.vcf.gz \
        $BASEDIR/hg19/chr$chr.dose.vcf.gz
    if [ $? -eq 0 ]
    then
        echo "Finished modifying header of chromosome $chr"
    else
        echo "Could not modify header of chromosome $chr"
        exit 1
    fi
    echo "Removing unsequenced samples from chromosome $chr"
    # the pysam script that I'm using appears to require a bcf file
    # all subsequent files will be gzipped vcf files, because they take up about half of the space
    # they can also be indexed with tabix, tho I'm using bcftools for this
    bcftools view -Ob \
        -s `cut -f1 ~/LabNotes/VCFindex.txt | cut -d/ -f 1 | cut -d- -f 1 | sort | uniq | join -t'|' - ~/LabNotes/mappings.txt | paste -s -d,` \
        -o $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.vcf.gz \
        $BASEDIR/hg19/chr$chr.dose.rename.vcf.gz
    if [ $? -eq 0 ]
    then
        echo "Finished removing unsequenced samples from chromosome $chr"
        rm header_temp.txt
    else
        echo "Could not remove unsequenced samples from chromosome $chr"
        rm header_temp.txt
        exit 1
    fi
    echo "Filtering SNPs with alt allele probabilities < 0.9 for $chr"
    python ~/LabNotes/Python/FilterVCF.py $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.vcf.gz \
      | bcftools view -Oz -o $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.filter_sites.vcf.gz
    if [ $? -eq 0 ]
    then
        echo "Finished filtering SNPs for chromosome $chr"
    else
        echo "Could not filter SNPs for chromosome $chr"
        exit 1
    fi
    echo "Adding rsIDs to $chr"
    bcftools index /c8000xd3/rnaseq-heath/Genotypes/Imputation3/hg19/chr$chr.dose.rename.filter_samples.filter_sites.vcf.gz
    bcftools annotate -c ID -Oz \
        -a /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/hg19/All_20161121.vcf.gz \
        -o $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.filter_sites.rsID.vcf.gz \
        $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.filter_sites.vcf.gz
    bcftools index $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.filter_sites.rsID.vcf.gz
    if [ $? -eq 0 ]
    then
        echo "Finished adding rsIDs to chromosome $chr"
    else
        echo "Could not add rsIDs to chromosome $chr"
        exit 1
    fi
    echo "Liftover of chromosome $chr from hg19 to GRCh38"
    bash ~/LabNotes/SubmissionScripts/Liftover.sh $chr
    if [ $? -eq 0 ]
    then
        echo "Finished liftover for chromosome $chr"
    else
        echo "Could not do liftover for chromosome $chr"
        exit 1
    fi
    python ~/LabNotes/Python/FilterVCF.py $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.vcf.gz \
      | bcftools view -Ob -o $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.filter_sites.vcf.gz
    rm header_temp.txt
done