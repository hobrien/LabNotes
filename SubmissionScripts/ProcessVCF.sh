#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

chr=$1
BASEDIR=/c8000xd3/rnaseq-heath/Genotypes/Imputation3

echo "Starting to process VCF for chr$chr"
ls /c8000xd3/rnaseq-heath/Mappings | cut -d- -f 1 | sort | uniq > ~/LabNotes/mappings.txt

if [ ! $BASEDIR/hg19/chr$chr.dose.rename.vcf.gz ]
then 
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
    rm header_temp.txt
fi

if [ ! $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.vcf.gz ]
then
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
fi

if [ ! $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.filter_sites.vcf.gz ]
then
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
fi

if [ ! $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.filter_sites.rsID.vcf.gz ]
then
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
fi

if [ ! $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.vcf.gz ]
then
    echo "Adding the stupid 'chr' to chromosome numbers for $chr"
    bcftools view $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.filter_sites.rsID.vcf.gz \
      | perl -pe 's/^(##contig=<ID=)?(\d+)/$1chr$2/' \
      | bgzip \
      > $BASEDIR/hg19/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.vcf.gz
    if [ $? -eq 0 ]
    then
        echo "Finished adding the stupid 'chr' to chromosome numbers for $chr"
    else
        echo "Could add the stupid 'chr' to chromosome numbers for $chr"
        exit 1
    fi
fi

if [ ! $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.vcf.gz ]
then
    echo "Liftover of chromosome $chr from hg19 to GRCh38"
    bash ~/LabNotes/SubmissionScripts/Liftover.sh $chr
    if [ $? -eq 0 ]
    then
        echo "Finished liftover for chromosome $chr"
    else
        echo "Could not do liftover for chromosome $chr"
        exit 1
    fi
fi

if [ ! $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.gz ]
then
    echo "Sorting GRCh38 VCF for $chr"
    bash ~/LabNotes/SubmissionScripts/SortVCF.sh /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.vcf.gz
    if [ $? -eq 0 ]
    then
        echo "Finished sorting GRCh38 VCF for $chr"
    else
        echo "Could not sort GRCh38 VCF for $chr"
        exit 1
    fi
fi

if [ ! $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.check.ref ] | \
   [ ! $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.check.nonSnp ] | \
   [ ! $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.check.dup ]
then
    echo "Running checkVCF on GRCh38 VCF for $chr"
    bash ~/LabNotes/SubmissionScripts/checkVCF.sh /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.gz
    if [ $? -eq 0 ]
    then
        echo "Finished running checkVCF on GRCh38 VCF for $chr"
    else
        echo "Could not run checkVCF on GRCh38 VCF for $chr"
        exit 1
    fi
fi
 
if [ ! $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.vcf.gz ]
then
    echo "Removing non-reference bases from GRCh38 VCF for $chr"
    cut -f 2 $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.check.nonSnp > $BASEDIR/GRCh38/chr$chr.excludedSNPs.txt
    cut -f 2 $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.check.ref | cut -f 2 -d: >> $BASEDIR/GRCh38/chr$chr.excludedSNPs.txt
    exset=`sort $BASEDIR/GRCh38/chr$chr.excludedSNPs.txt |uniq | perl -pe 's/^/POS=/' | paste -s -d'|'`
    if [[ ${#exset} > 0 ]]
    then
        bcftools filter -Oz -e $exset $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.gz \
          > $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.vcf.gz
    else
        cp $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.gz \
           $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.vcf.gz
    fi
    if [ $? -eq 0 ]
    then
        echo "Finished removing non-reference bases from GRCh38 VCF for $chr"
    else
        echo "Could not remove non-reference bases from GRCh38 VCF for $chr"
        exit 1
    fi
fi

if [ ! $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.filter_dup.vcf.gz ]
then
    echo "Removing duplicated positions from GRCh38 VCF for $chr"
    cut -f 2 $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.check.dup | cut -f 2 -d: >> $BASEDIR/GRCh38/chr$chr.excludedSNPs.txt
    exset=`sort $BASEDIR/GRCh38/chr$chr.excludedSNPs.txt |uniq | perl -pe 's/^/POS=/' | paste -s -d'|'`
    if [[ ${#exset} > 0 ]]
    then
        bcftools filter -Oz -e $exset $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.vcf.gz \
          > $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.filter_dup.vcf.gz
    else
        cp $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.vcf.gz \
           $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.filter_dup.vcf.gz
    fi
    if [ $? -eq 0 ]
    then
        echo "Finished removing duplicated positions from GRCh38 VCF for $chr"
    else
        echo "Could not remove duplicated from GRCh38 VCF for $chr"
        exit 1
    fi
fi


if [ ! $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.filter_dup.vcf.check.ref ] | \
   [ ! $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.filter_dup.vcf.check.nonSnp ] | \
   [ ! $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.filter_dup.vcf.check.dup ]
then
    echo "Rerunning checkVCF on filtered GRCh38 VCF for $chr"
    bash ~/LabNotes/SubmissionScripts/checkVCF.sh /c8000xd3/rnaseq-heath/Genotypes/Imputation3/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.filter_dup.vcf.gz
    if [ $? -eq 0 ]
    then
        echo "Finished rerunning checkVCF on filtered GRCh38 VCF for $chr"
    else
        echo "Could not rerun checkVCF on filtered GRCh38 VCF for $chr"
        exit 1
    fi
fi

if [ ! $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.filter_dup.vcf.gz.csi ] 
then
    echo "Indexing (csi) filtered GRCh38 VCF for $chr"
    bcftools index $BASEDIR//GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.filter_dup.vcf.gz
    if [ $? -eq 0 ]
    then
        echo "Finished indexing (csi) filtered GRCh38 VCF for $chr"
    else
        echo "Could index (csi) filtered GRCh38 VCF for $chr"
        exit 1
    fi
fi

if  [ ! $BASEDIR/GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.filter_dup.vcf.gz.tbi ]
then
    echo "Indexing (tbi) filtered GRCh38 VCF for $chr"
    tabix $BASEDIR//GRCh38/chr$chr.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.filter_nonSNP.filter_dup.vcf.gz
    if [ $? -eq 0 ]
    then
        echo "Finished indexing (tbi) filtered GRCh38 VCF for $chr"
    else
        echo "Could index (tbi) filtered GRCh38 VCF for $chr"
        exit 1
    fi
fi

echo "Finished to processing VCF for chr$chr"

