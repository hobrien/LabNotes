
Table of Contents
=================

  * [Table of Contents](#table-of-contents)
  * [Create VCF with imputed SNPs](#create-vcf-with-imputed-snps)
    * [Prepare files for imputation](#prepare-files-for-imputation)
    * [Concatenate and filter imputed data](#concatenate-and-filter-imputed-data)
      * [bcftools](#bcftools)
    * [Add SNP IDs to imputed VCFs](#add-snp-ids-to-imputed-vcfs)
  * [Prepare Sequencing data](#prepare-sequencing-data)
    * [Sequence data QC](#sequence-data-qc)
      * [FastQC](#fastqc)
    * [Trimming sequences](#trimming-sequences)
      * [Trim Adaptors from reads using cutadapt](#trim-adaptors-from-reads-using-cutadapt)
      * [Use trim\_galore cutadapt wrapper](#use-trim_galore-cutadapt-wrapper)
  * [Mapping](#mapping)
    * [Start mapping reads with tophat](#start-mapping-reads-with-tophat)
    * [Switch to mapping with HISAT](#switch-to-mapping-with-hisat)
    * [Mapping QC](#mapping-qc)
    * [Try to run RNAseq\-specific QC:](#try-to-run-rnaseq-specific-qc)
      * [RNA\-SeQ](#rna-seq)
      * [CollectRnaSeqMetrics](#collectrnaseqmetrics)
      * [RSeQC](#rseqc)
    * [I am going to need to do a liftover before running this because the vcf is for hg19](#i-am-going-to-need-to-do-a-liftover-before-running-this-because-the-vcf-is-for-hg19)
      * [LiftoverVcf (picard\-tools)](#liftovervcf-picard-tools)
      * [CrossMap](#crossmap)
    * [WASP](#wasp)
    * [Clip overlapping reads](#clip-overlapping-reads)
  * [eQTL analyses](#eqtl-analyses)
    * [MatrixEQTL](#matrixeqtl)
  * [Transcript Identification](#transcript-identification)
    * [Cufflinks](#cufflinks)
    * [Cuffcompare](#cuffcompare)
    * [Make DB of GencodeGTF:](#make-db-of-gencodegtf)
    * [Run Cufflinks on SRA data from Liu:2016ji (GSE71315)](#run-cufflinks-on-sra-data-from-liu2016ji-gse71315)
  * [Expression analysis](#expression-analysis)
  * [Cell\-type deconvolution](#cell-type-deconvolution)
  * [Differential Expression Analyses](#differential-expression-analyses)
    * [HTSeq\-count](#htseq-count)
    * [DESeq2/EdgeR](#deseq2edger)
    * [DEXSeq](#dexseq)
    * [JunctionSeq](#junctionseq)
    * [<a href="https://pachterlab\.github\.io/kallisto/">Kallisto</a>](#kallisto)
    * [Cufflinks](#cufflinks-1)
  * [SNP calling](#snp-calling)
  * [<a href="https://software\.broadinstitute\.org/gatk/gatkdocs/org\_broadinstitute\_gatk\_tools\_walkers\_rnaseq\_ASEReadCounter\.php">ASEReadCounter</a>](#asereadcounter)

Created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc.go)

# Create VCF with imputed SNPs
## Prepare files for imputation
- All steps needed for this analysis are run by the check-bim script on [this](http://www.well.ox.ac.uk/~wrayner/tools/) site (v4.2.5, with a few modifications for additional processing steps)
    - changed output from bed to VCF
    - added sort and compression steps
    - final output written to subfolder, which makes upload to [Michigan imputation server](https://imputationserver.sph.umich.edu/start.html) easier
- Analysis steps (this has to be done on the iMac because it requires ca. 12 GB of RAM):
    - ```cd BTSync/FetalRNAseq/Genome-wide\ genotyping/```
    - ```plink --bfile FB_Merged --freq```
    - ```mkdir Imputation3```
    - ```cd Imputation3```
    - ```perl ../LabNotes/Perl/HRC-1000G-check-bim.pl -b ../FB_Merged.bim -f ../plink.frq -r ../../Ref/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h```
    - ```bash Run-plink.sh```
    - I've now written a script to run all of these steps on the server (ProcessPLINK.sh)
- Run imputation using HRC r1.1 2016 and Eagle v2.3 for phasing with Mixed population for QC (see ScreenShots/Imputation3.png)
    - files are uploaded from BTSync/FetalRNAseq/Genome-wide\ genotyping/Imputation3
    - Imputed data are in /Volumes/FetalRNAseq/ImputedGenotypes/Imputation3 and in /c8000xd3/rnaseq-heath/Genotypes/Imputation3 on rocks
        
## Concatenate and filter imputed data
### bcftools
- modify vcf headers to add info about all filter classes (PASS/Genotyped/ Genotyped_only) (see [this](https://github.com/samtools/bcftools/issues/470) post)
    - ```head -13 ~/LabNotes/header.txt > header_temp.txt```
    - ```echo "##contig=<ID=$chr>" >> header_temp.txt```
    - ```tail -1 ~/LabNotes/header.txt >> header_temp.txt```
    - ```bcftools reheader -h header_temp.txt -o /c8000xd3/rnaseq-heath/Genotypes/Imputation3/hg19/chr$chr.dose.rename.vcf.gz /c8000xd3/rnaseq-heath/Genotypes/Imputation3/hg19/chr$chr.dose.vcf.gz```
- Change the names of the samples in the VCF headers to match the samples (this step is unnecessary if I use the new names when I reheader in the step above)
    - ```cut -f 1 ~/LabNotes/VCFindex.txt | cut -d'/' -f 1 | cut -d'-' -f 1 > ~/LabNotes/VCFindex2.txt```
    - ```for chr in {1..22}; do bcftools reheader -s ~/LabNotes/VCFindex2.txt -o /c8000xd3/rnaseq-heath/Genotypes/Imputation3/hg19/chr$chr.dose.rename.vcf.gz /c8000xd3/rnaseq-heath/Genotypes/Imputation3/hg19/chr$chr.dose.vcf.gz; done```
- Filter out samples that we do not have sequence data for
    - ```ls /c8000xd3/rnaseq-heath/Mappings > ~/LabNotes/mappings.txt```
    - ```for chr in {1..22}; do bcftools view -s `cut -f1 ~/LabNotes/VCFindex.txt | cut -d/ -f 1 | cut -d- -f 1 | sort | join -t'|' - ~/LabNotes/mappings.txt | paste -s -d,` -o  /c8000xd3/rnaseq-heath/Genotypes/Imputation3/hg19/chr$chr.dose.rename.filter_samples.vcf.gz /c8000xd3/rnaseq-heath/Genotypes/Imputation3/hg19/chr$chr.dose.rename.vcf.gz; done```

- I need to combine into a single file
    - ```find Raw_output -name chr\*.dose.vcf.gz | awk '{ print length, $0 }' |sort -n -s | cut -d" " -f2- > vcf_files.txt```
    
    - ```bcftools concat -o All_chromosomes.vcf.gz -f vcf_files.txt -O b```
- I also need to pull out data on 2 individuals for Nick's collaborator
- Added database table (PC_analysis) where Sentrix_Full matches VCF IDs and BrainBankID matches 
    - ```echo {1..22} |xargs -n 1 -I % bcftools view -s 17_9702504079_R09C01,67_9702504147_R10C01 -O b -o Subset/Chr_%.vcf Raw_output/chr_%/chr%.dose.vcf.gz```
    - ```ls -rth Subset/ |perl -pe 's/^/Subset\//' >subset.txt```
        - this sorts correctly because it's the order in which they were analysed. Find puts 3 after 33
    - ```bcftools concat -O b -f subset.txt -o subset.vcf.gz ``` 
    - ```python ~/BTsync/FetalRNAseq/LabNotes/Python/ChangeSampleID.py subset.vcf.gz```
    - ```rm temp.head```
- replace genotyping IDs with Brain Bank IDs
    - ```python ~/BTsync/FetalRNAseq/LabNotes/Python/ChangeSampleID.py All_chromosomes.vcf.gz```
    - ```bcftools reheader -h temp.head -o All_chromosomes_relabelled.vcf.gz All_chromosomes.vcf.gz```
    - ```rm temp.head```

## Add SNP IDs to imputed VCFs
- SNP DB [Schema](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp146.sql) and [data](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp146.txt.gz) downloaded from the [USCSC Genome Browser](https://genome.ucsc.edu/) and imported into FetalRNAseq mySQL DB
    - ```mysql -u root FetalRNAseq < SQL/snp146.sql```
    - chromosome (field 1), position (2) and ref (4) and alt (5) alleles can be used to uniquely identify each SNP (I hope!)
        - there can be several SNPs with the same start position, but they appear to have different end positions
        - there can also be > 2 observed variants at a position                                
        - ```grep -v '^#' FB_Merged_chr22.vcf | cut -f 1,2,3 |python ../LabNotes/Python/AddSNPID.py | wc -l```
        - all rsIDs match, except for some that have been merged together
        - there are also a lot of duplicates, that presumably will be merged in the future
        - if I look up SNPs by start position and end position, it should be posible to get a single SNP (or multiple duplicate SNPs) (I think!).
        - comparing the observed alleles in the DB to the ref/alt alleles in the VCF is a pain because
            - a) the observed alleles can be on the opposite strand to the ref Allele
            - b) there can be more than 2 observed alleles
        - there is a exceptions field that will indicate if there are duplicates.
        - All of the warnings that I get on chromosome 22 are now due to SNPs that have been merged
    - I've now modified AddSNPID.py to do what the name suggests:
        - ```bcftools view All_chromosomes_relabelled.vcf.gz | python ../LabNotes/Python/AddSNPID.py |bgzip -c >All_chromosomes_rsID.vcf.gz```
        
# Prepare Sequencing data
## Sequence data QC
### FastQC        
- Running FastQC on Edinburgh data
    - I would like to 
        - (a) work with compressed data and 
        - (b) pass the file name to the submission script so i can do something like ``` find . fastq.gz | xargs -n 1 qsub SubmissionScripts/FastQC.sh```
    - Neither of these is working ATM
        - the [recommended](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt) method of working with compressed data (```zcat seqs.fastq.gz | fastqc stdin```) produces a file called stdin_fastqc.zip which can't be opened
        - it appears to be possible to just supply the name of the compressed file though
        - $1 and $@ don't work for the input because they get split on spaces in the filenames
        - '$@' isn't recognised, but "$@" seems to work great
    - ```find /c8000xd3/databank/foetal-rna/1st\ batch\ Edinburgh\ Sequencing/ -name *.sanfastq.gz -print0 |xargs -0 -n 1 qsub SubmissionScripts/FastQC.sh```
    - ```find /c8000xd3/databank/foetal-rna/Edinburgh\ 2nd\ batch\ sequencing/ -name *.sanfastq.gz -print0 |xargs -0 -n 1 qsub SubmissionScripts/FastQC.sh```
    - ```find /c8000xd3/databank/foetal-rna/Exeter\ sequencing/ -name \*fastqc\* -print0 | xargs -n 1 bash SubmissionScripts/cp.sh```
    - FastQC.sh:
        - ```~/src/FastQC/fastqc --outdir=/home/heath/FastQC "$@"```
    - cp.sh:
        - ```cp "$@" >FastQC/```
    - Uncompress all files and concatenate
        - ``` unzip -d FastQC/Uncompressed 'FastQC/*.zip'```
        
        - ``` find Uncompressed/ -name summary.txt | grep -v SRR | grep -v Undetermined | grep -v _CER_ | grep -v _T_ctx_ | grep -v trim | xargs perl -pe 's/.*\///' >summary.txt```
    - Extract the number of reads for each SampleID
        - ```find Uncompressed/ -name fastqc_data.txt | xargs grep 'Total Sequences' | grep -v trim | grep -v SRR | grep -v Undetermined | grep -v _CER_ | grep -v _T_ctx_ | perl -pe 's/(\.sanfastq)?_fastqc\/fastqc_data\.txt\:Total Sequences//' | perl -pe 's/.*\///' > seq_lengths.txt```
    - Results are analysed in FastQC.md

## Trimming sequences
### Trim Adaptors from reads using cutadapt
- I was able to install cutadapt on rocks by first installing pip in ~/.local (``` python ~/src/get-pip.py --user```) using it:
    - ```~/.local/bin/pip install --user --upgrade cutadapt```
- After much struggle with the stupid space names, I finally had the genius of symlinking the data in my home folder:
    - find /c8000xd3/databank/foetal-rna/ -name \*fastq.gz |python ~/SubmissionScripts/LinkRaw.py
- Annoyingly, must of the Exeter data is uncompressed, and now can't (easily) be changed. I will just compress them and store them along with the symlinks in my home folder
    - ```find /c8000xd3/databank/foetal-rna/Exeter\ sequencing/ -name *.fastq | grep -v 17921-l1_CGATGT_L006_R1_001 | python SubmissionScripts/LinkRaw.py```
        - grep command is to skip one file that I tried by itself
    - ```find ~/Temp/ -name *.fastq | grep -v 17921-l1_CGATGT_L006_R1_001 | xargs -n 4 qsub ../SubmissionScripts/Compress.sh```
    - ```rm -r ~/Temp```
- Make sure all files are matched
    - ```find ../Raw/ -name *.fastq.gz | sort | xargs -n 2 python ../SubmissionScripts/CheckNames.py```
        - no warnings
    - ```find ../Raw/ -name *.fastq.gz | sort | xargs -n 2 | wc -l```
        - 81 (19 Exeter + 32 Edinburgh1 + 30 Edinburgh2)
- Trim Edinburgh data while Exeter data is still compressing
    - Edinburgh data only: ```../Raw/ -name \*TP\*.fastq.gz | sort``
    - Exeter data only ```find ../Raw/ -name \*TP\*.fastq.gz | sort```
    - ```find ../Raw/ -name \*TP\*.fastq.gz | grep -v 150429_D00200_0258_BC6UARANXX_4_IL-TP-0 | sort |xargs -n 2 qsub ../SubmissionScripts/CutAdapt.sh```
        - this skips the 2 that are already trimming
- Rerun FastQC on trimmed data
    - ```find ../Trimmed -name *.fastq.gz -print0 |xargs -0 -n 1 qsub ../SubmissionScripts/FastQC.sh```
- download and analyse:
    - ```unzip -d FastQC/Uncompressed 'FastQC/*_trimmed_fastqc.zip'```
    - ```find Uncompressed/ -name summary.txt |grep trimmed |xargs perl -pe 's/_trimmed.*//' >>trimmed_summary.txt```

### Use trim_galore cutadapt wrapper
- For some obscure reason, the reads are becoming unpaired after trimming. I get very high pair congruence when mapping the raw data with Tophat, but very low congruence with the trimmed data
- I've also decided that it's a good idea to do some quality trimming because unlike Star, Tophat does not do read clipping
- [trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore) is a perl wrapper for cutadapt that makes it a little more flexible and easy to use. It also solves this problem
    - After doing the mapping with quality trimming and adapter trimming, everything looks good except that there's a big spike in SNP frequency at the 5' end of the reads.
    - I am going to solve this by clipping 5 bp from the 5' end of each read in addition to the quality and adapter trimming:
        - ```find /c8000xd3/databank/foetal-rna/Exeter_sequencing/15533_Oct_2014/ -name *.fastq.gz | xargs qsub SubmissionScripts/Trim.sh```

# Mapping            
## Start mapping reads with tophat
- follow the instructions [here](http://www.illumina.com/documents/products/technotes/RNASeqAnalysisTopHat.pdf) to get started with tophat
    - ```wget --ftp-user=igenome --ftp-password=G3nom3s4u ftp://ftp.illumina.com/Homo_sapiens/NCBI/GRCh38Decoy/Homo_sapiens_NCBI_GRCh38Decoy.tar.gz```
    - ```tar -xzf Homo_sapiens_UCSC_hg19.tar.gz```
    - ```tophat --keep-fasta-order --transcriptome-index /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.inx --library-type fr-secondstrand --mate-inner-dist 150  --mate-std-dev 50 --num-threads 8 --output-dir /home/heath/Mappings/15533_150 /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome /home/heath/Trimmed/15533_TGACCA_L007_R1_001_trimmed.fastq.gz /home/heath/Trimmed/15533_TGACCA_L007_R2_001_trimmed.fastq.gz```
## Switch to mapping with HISAT
- this gives very similar results, but runs a lot faster
- I'm using Tophat mapping for the differential expression analysis, but HISAT for Allele-specific analyses

## Mapping QC
- install [bamQC](https://github.com/s-andrews/BamQC) on rocks
    -``` cd src```
    - ```git clone https://github.com/s-andrews/BamQC```
    - ```cd BamQC```
    - ```ant```
    - ```chmod 700 bin/bamqc```
    - ```ln -s ~/src/BamQC/bin/bamqc ~/bin/bamqc```
        
- [download](http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format) a complete GTF file for hg19:
    - ```curl http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf >~/bin/genePredToGtf```
    - ```chmod +x ~/bin/genePredToGtf```
    - create file ~/.hg.conf:
        ```
        db.host=genome-mysql.cse.ucsc.edu
        db.user=genomep
        db.password=password
        central.db=hgcentral
        ```
    - ``` chmod 600 ~/.hg.conf```
    - ```genePredToGtf hg19 knownGene ~/Index/knownGene.gtf```
- run bamQC
    - ```bamqc --gff ~/BTSync/FetalRNAseq/Reference/knownGene.gtf 15533/accepted_hits.bam```
- plot insert sizes
    - ```cat BamQC/15533_300/accepted_hits_bamqc/bamqc_data.txt | python LabNotes/Python/ExtractInsertSize.py >BamQC/15533_300/insert_sizes.txt```    
    - results are analysed in BamQC.md

## Try to run RNAseq-specific QC:
### RNA-SeQ
- Eilis used something called RNA-SeQ. She went thru several steps to format the mapping before this would work 
    - ```java -jar ~/Documents/src/picard-tools-1.119/CreateSequenceDictionary.jar R=genome.fa O=genome.bam```
    - ```java -jar ~/Documents/src/picard-tools-1.119/AddOrReplaceReadGroups.jar I=accepted_hits.bam O=accepted_hits2.bam RGLB="totalRNA" RGPL="Illumina" RGPU="1" RGSM="ID_15533"```
    - ``` samtools index accepted_hits2.bam```
    - ```java -jar ~/Documents/src/picard-tools-1.119/ReorderSam.jar I=accepted_hits2.bam O=accepted_hits3.bam R=~/BTSync/FetalRNAseq/Reference/genome.fa```
    - ```samtools index accepted_hits3.bam```
    - ```java -jar ~/Documents/src/RNA-SeQC_v1.1.8.jar -n 1000 -s "ID_15533|accepted_hits3.bam|Test" -t ~/BTSync/FetalRNAseq/Reference/gencode.v19.annotation.gtf -r ~/BTSync/FetalRNAseq/Reference/genome.fa -o RNAseQC -gc ~/BTSync/FetalRNAseq/Reference/gencode.v7.gc.txt```
        - still not working 
        
### CollectRnaSeqMetrics
- Run [CollectRnaSeqMetrics]( http://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics) from [Picard Tools](http://broadinstitute.github.io/picard/)
            
    - had to install a more recent version of java on rocks
            
    - downloaded a correctly formatted [file](https://gist.github.com/slowkow/b11c28796508f03cdf4b) with rRNA intervals     
        - I'm pretty sure this stat is not correct. I don't know if it's a problem with this file or with the command I used
    - I had to build a dictionary for the reference Seq
        - ```/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=/home/heath/Ref/hg19.fa O=/home/heath/Ref/hg19.dict```    
    - A (sort of) explanation of the output is [here](https://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics)

### RSeQC
- Run [RSeQC](http://rseqc.sourceforge.net)
    - I figured it would be best to use a BED file made from the actual annotation I'm using
        - [BEDOPS](https://bedops.readthedocs.org/en/latest/) seemed like a good option for this, but I can't figure out how to get a proper [BED12](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file from it, so I grabbed the recommended [perl script](https://code.google.com/p/ea-utils/source/browse/trunk/clipper/gtf2bed)
    - These scripts are working great, but I'm getting very few validly paired reads, even when using the fr-secondstrand option.
        - This affects the stats, and also possibly the number of pairs used for the inner_distance info (only 1800 of 1,000,000)
        - I'm not running tophat in all possible strand modes to pin down this issue. I suspect the problem is the insert size though                
    - Get Chromosome sizes
        - ```bash Lab_notes/Bash/fetchChromSizes hg19 >Reference/hg19.chrom.sizes```
    - ```bam_stat.py -i ~/Documents/Mappings/15533_300/accepted_hits.bam```
    - ```inner_distance.py -i accepted_hits.bam -o ~/BTSync/FetalRNAseq/BamQC/15533_1000/15533_1000 -r ~/BTSync/FetalRNAseq/Reference/hg19_RefSeq.bed -u 5000 -s 50 > /dev/null```    
    - Check for mapping to rRNA
        - It looks like the ribosomal RNA genes are not annotated in GRCh38, but if I blast the repeat unit (       /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/AbundantSequences/humRibosomal.fa) against the genome, I get a series of hits with scores > 7,000 in the following regions:
            - chr21:8250197-8472360 (222 kb)
            - chr22_KI270733v1_random:122273-179772 (57 kb)
            - chrUn_GL000220v1:105424-161802 (56 kb)
                - ```cat humRibosomal.bl |awk '{ if ($12 >= 7000) print $0 }'```
        - Convert blast results to BED, then run split_bam.py:
            - ```cat /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/AbundantSequences/humRibosomal.bl |blast2bed12.py > /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/AbundantSequences/humRibosomal.bed````
    - This is working great now, except that insertion_profile.py was reporting results for 52 bp reads. I'm testing it now with the -l 100 option to see if that helps
    - It would be nice if it could separate SNPs from reads 1 and 2, but otherwise, it's pretty nice                     
- for HISAT mappings, I've not bothered to run the full suite of QC metric, but I am calculating basic stats because I find it to be a lot more useful than flagstat
    - it reports the total number of records, but then breaks them down into unmapped, unique, non-unique, and secondary alignments
    - taking a sum of the first three categories gives the total number of input reads
    - when trimming, the total number of reads remains constant, but some reads end up with a length of zero
    - the number of reads reported in this summary is lower than the number of reads in the input file (by 4.4 million), so the mapper must be trowing away short reads
    - 3.7 million less reads are unmapped after trimming, and 2.9 million less are unique while 2.2 million more are non-unique
    - this means basically means that 3.7 crap reads that don't map plus 700 k crap reads that do map are filtered out and that 2.2 million uniquely mapping reads are non-unique after trimming
    - this means that the difference in the number of records (1.9 million) is less than the difference in the number of reads (4.4 million), because 2.2 million non-unique reads map a total of 4.7 million times
    - while 2.9 million fewer reads are uniquely mapped, 3.5 million fewer reads are mapped in proper pairs. At first I thought this mean that the mapping accuracy is higher with the untrimmed data, but each read that is unmapped on non-uniquely mapped results in two reads being unpaired

## I am going to need to do a liftover before running this because the vcf is for hg19
### LiftoverVcf (picard-tools)    
- there is a picard tool for this purpose. I just need to download the chain from UCSC
    - ```~/bin/java -jar ~/src/picard-tools-2.1.1/picard.jar LiftoverVcf```
        - this isn't working. all SNPs are being rejected. I think this is because the chromosomes are coded '1', etc rather than 'chr1', etc.
        - I am attempting to fix this:
            - ```bcftools view chr1.dose.vcf.gz |perl -pe 's/^1/chr1/' | bgzip >chr1.recoded.vcf.gz```
                - I had to install a newer version of [bcftools](http://www.htslib.org/download) on the server because the command line interface has changed
    - this is a right pain. I'm having major issues with b37 (lacks 'chr') and hg19 (has chr). See [this](http://gatkforums.broadinstitute.org/gatk/discussion/63/errors-about-input-files-having-missing-or-incompatible-contigs).
    - the UCSC chain file (hg19ToHg38.over.chain) has the 'chr' prefix while the Ensembl one (GRCh37_to_GRCh38.chain) works with or without the prefix, but outputs contig names without them.
    - I recoded the fasta file (genome_recoded.fa) to work with the Ensembl chain file, but it is going to cause downstream problems because the BAM files have the prefixes
    - there is an additional problem in that picard-tools-2.1.1 is unable to find the sequences even when the names match. something about [using bytes instead of names](https://github.com/broadinstitute/picard/releases/download/2.2.0/README.txt). Newer versions work
    - finally, the memory requirements for this are HUGE! It took 12 GB to process 20,000 lines. Chromosome 1 has 3.7 million lines, so that's going to require truly massive memory, or a hell of a lot of splitting and concatinating.

### CrossMap
- fortunately, it looks like something called [CrossMap](http://crossmap.sourceforge.net) is going to come to my rescue. 
    - It works fine on 100,000 lines. We'll see how it does on all of them
    - The only problem is that it produces an uncompressed VCF, so I'll have to run bgzip after it finishes    

## WASP
- Run [WASP](https://github.com/bmvdgeijn/WASP) allele-aware alignment
    - I compiled snp2hd5, which worked fine, though there was a fairly obvious problem with the vcf parsing that I had to be fixed
    - It is reporting an error on the last line of the chr22 vcf
        - This is an offset error that I think it related to the size of the matrix, tho the code looks fine to me
        - No idea if this is going to cause problems or not (beyond meaning that this snp in not included)
        - It is also reporting that the data are unphased, which appears to be correct. No idea how the imputation could have happened without phasing the data though.
        - VCF files are improperly sorted for some reason:
            - ```/share/apps/vcftools-0.1.14/bin/vcf-sort chr2.GRCh38.vcf.gz | bgzip > chr2.GRCh38.sort.vcf.gz```
            - That would be because I selected the unphased output option
    - Need a list of SNPs for each chromosome for the mapping:
        - ```bcftools view -H chr1.GRCh38.vcf.gz cut -f1,2,4,5 | gzip >SNPs/chr1.snps.txt.gz```
    - This gives a list of all SNPs in all samples. It might make more sense to use sample-specific lists of SNPs: 
        - ```bcftools view -H chr1.GRCh38.vcf.gz | cut -f1,2,3,4,5,6,7,8,9,10 | grep -v '0|0' | cut -f1,2,4,5 | gzip > 15240/chr1.snps.txt```
    - This gives 309k SNPs as compared to 3.7M  
    - Get index position for each sample in VCF files:
        - ```bcftools query -l ../../ImputedGenotypes/Imputation2/chr_1/chr1.dose.vcf.gz | head -96 | python GetVCFindex.py > ../VCFindex.txt```  
    - This tool has been [updated](https://github.com/bmvdgeijn/WASP/releases/tag/v0.2), including a bunch of fixes to things that have been vexing me
        - uses HDF5 to store SNP info instead of hokey text file
        - considers only SNPs that are polymorphic in the sample
        - "fixed several mapping pipeline bugs related to paired-end reads" (hopefully including the problems with mates being dropped and reads being duplicated)
        - some other fixes that will hopefully solve the taking forever problem
        - I had to update a number of things in anaconda to get this to work properly
            - ```conda install -c anaconda gmp=5.1.2```
            - ```conda update pysam```
        - Working on generating HDF5 from VCF, but I'm getting segfaults for a couple of the VCF files (/c8000xd3/rnaseq-heath/Genotypes/Imputation2/Sorted/chr14.GRCh38.sort.vcf.gz, /c8000xd3/rnaseq-heath/Genotypes/Imputation2/Sorted/chr19.GRCh38.sort.vcf.gz)
            - There is also a problem with chr22 that gives an error saying "SNP position (89346516) is outside of chromomosome  range:1-72058697861300481"
            - I have no idea why it gives the range as 1-72058697861300481, which is 25 million x the total size of the genome, but 89346516 is outside the range of chr22 (50818468)
            - The VCF for chr22 actually goes up to position 154768956, which is more than the length of all be ut the first 7 chromosomes
            - the name of the SNP at this position is chr22:18047467, which refers to the position in the hg19 build
            - rs188965487 is at this position in hg19, but the position is listed as "no hit" in GRCh38.p2
            - rs575160859 is the last SNP in chr22.GRCh38.vcf.gz (unsorted). It's position is listed as 50805809, which is it's position in GRCh38.p2. It's name is listed as 22:51244237, which corresponds to it's position in hg19
            - It very much looks like most of the SNPs in chr22.GRCh38.vcf.gz are correct, and in the correct order, but for some reason, some SNPs that are not present in GRCh38 were given extremely high coordinates that are throwing up sort order issues and are outside of the range of the chromosome
            - There are a total of 8 SNPs that have coordinates beyond the length of the chromosome, 5 within 44 bp in hg19 and 3 others within 241 bp. I haven't checked all of them, but both of these clusters have SNPS with no hit in GRCh38
            - These are easy enough to exclude, but what I really need is a sanity check that ensures that all SNPs in the lift-over file have a matching reference base at the listed position in GRCh38
            - FixVCF.py flags 7 of these as 'FailedGetBase' in file called test.check.ref (absolutely no idea why it missed chr22:89346549)
            - There is also a SNP in the same file flagged as MismatchRefBase, but it's a mismatch in both the original VCF and the liftover. Looks like a case of a strand swap during imputation
            - It also flags 4000 duplicate sites, which are different alternate alleles at the same position
            - There are also complaints about Alternative alleles with freq > 50% and monomorphic sites, but I think I can safely ignore them
            - Filtering on position appears to clean everything up:
                - ```bcftools filter -Oz -e 'POS>=89346516' chr22.GRCh38.vcf.gz > chr22.GRCh38.filter.vcf.gz```
                 
## Clip overlapping reads
- the tool of choice for this appears to be [clipOverlap](http://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap) from [bamUtil](http://genome.sph.umich.edu/wiki/BamUtil)
- this clips the read with the lowest quality score, which isn't as good as comparing the bases of the overlapping reads, but is probably good enough
- this turned out to be a huge pain to install, until [bioconda](https://bioconda.github.io/index.html#setup) saved the day
    - ```conda config --add channels r```
    - ```conda config --add channels bioconda```
    - ```conda install bamutil```
  
# eQTL analyses
## MatrixEQTL
- this requires 5 (!) files with bespoke formatting:
    - [expression file](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/GE.txt)
        - I can get normalised counts from the EdgeR analysis. I'll probably just code this in the R script as well
    - [genotypes file](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/SNP.txt)
        - I've written a script to convert a VCF file to this format (it prints SNP position file at the same time):
            ```bcftools view genotypes.vcf.gz | python ~/LabNotes/Python/VCF2numeric.py genotypes.txt snp_pos.txt``` 
    - [covariates file](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/Covariates.txt)
        - I added code to create this file from the target file and import it to the MatrixEQTL.R script
    - [gene location file](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/geneloc.txt)
        ```echo "geneid\tchr\ts1\ts2" > geneloc.txt```
        ```cat genes.gtf | awk '{if ($3 == "gene") print $10, $1, $4, $5}' | sed 's/[".;]//g' >> geneloc.txt```
    - [SNP location file](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/snpsloc.txt)
        - see genotypes file above
- all of this is now done in the MatrixEQTL.R script
                           
# Transcript Identification
## Cufflinks
- Run [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks)
    - Cufflinks is running on 15533, but it's taken 3 days so far and no indication of progress
    - Cufflinks finished after 1 day when I ran it on 15468, which isn't a much smaller dataset, so I'm not sure what's going on with that
    - Cufflinks includes an option (-M) to mask reads matching the specified features (ie; rDNA), but it's probably easier to just use the BAM files that already exclude these reads
    - I think it would also be helpful to get rid of all of the non-chromosome sequences.
        - FilterBAM.sh should do this

## Cuffcompare
- I downloaded a set of masked chromosome files from [UCSC](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz) which can be used to filter results with the -s option:
    - ```cuffcompare -V -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf -s /c8000xd3/rnaseq-heath/Ref/chroms```
- I've now run Cuffcompare on all samples
    - The input GTF files from Cufflinks have between 1.4 million and 2 million exons and between 300,000 and 1 million transcripts
        - ```ls | grep gtf$ | xargs -n 1 awk '{if ($3 == "exon") count ++} END{print count}' | sort -n```
        - ```ls | grep gtf$ | xargs -n 1 awk '{if ($3 == "transcript") count ++} END{print count}' | sort -n```
        - 17198 has 30% more exons and more than twice as many transcripts as any other sample
    - the combined GTF assembly 3.1 million exons (```wc -l Combined.combined.gtf```), 1.5 million transcripts (```wc -l Combined.tracking```) and 258,321 loci (```wc -l Combined.loci```)
- Concatinate all tmap results and upload them to DB
    - ```grep 0  *.tmap | perl -pe 's/Combined\.(\d+)\.gtf\.tmap:/$1\t/' > All.gtf.tmap```
- Get distribution of classes for all "multiple types" entries in Combined.tracking:
    -``` grep '\t\.\t'  Combined.tracking | python ../LabNotes/Python/MixedClasses.py >Mixed.txt```  
- A common novel transcript type is large single exon transcripts classed as 'i' or 'x'
    - these are very likely to represent unspliced introns and should probably be filtered out
    - does a pretty poor job of distinguishing x from i as far as I can tell
    - need to plot exon size by feature type to investigate this

## Make DB of GencodeGTF:
    - ```cat  gencode.v24.chr_patch_hapl_scaff.annotation.gtf |python ~/BTSync/FetalRNAseq/LabNotes/Python/GTF2CSV.py > ~/BTSync/FetalRNAseq/Reference/gencode.v24.chr_patch_hapl_scaff.annotation.csv```
    - ```mv features.csv ~/BTSync/FetalRNAseq/Reference/gencode.v24.chr_patch_hapl_scaff.features.csv```
    - ```mysql db_name < ~/BTSync/FetalRNAseq/LabNotes/SQL/refGTF.sql```
    - ```SELECT DISTINCT feature FROM GencodeFeatures```
        - transcript_support_level, exon_number, level, gene_type, transcript_type, gene_name, transcript_name, tag, ccdsid, exon_id, gene_id, protein_id, transcript_id, gene_status, transcript_status, havana_gene, havana_transcript, ont        

## Run Cufflinks on SRA data from Liu:2016ji (GSE71315)
- downloaded sra toolkit and ran fastq-dump (GetSRA.sh)
- hopefully I get similar results to what I have for our data
- I've downloaded all bulk tissue riboZero RNAseq datasets and Tophat is running on one
- The mapping QC looks VERY similar to our data, but a lot less reads, resulting in 10+ fold fewer novel junctions
    - This suggests that using a higher coverage threshold is the answer to cleaning up the novel cufflinks data
      
# Expression analysis
- Analyse expressed SNPs
    - Run mpileup on SNPs from grant:
        - ```cat ~/BTSync/FetalRNAseq/Info/ExpressedSNPs.txt | python ~/BTSync/FetalRNAseq/LabNotes/Python/GetSNPpos.py | xargs -n 1 -I % samtools mpileup -d 8000 -f ~/BTSync/FetalRNAseq/Reference/genome.fa -r % -ABQ 0 accepted_hits.bam |python ~/BTSync/FetalRNAseq/LabNotes/Python/CountBases.py ```  
- Concatinate tracking files and upload to DB:
    - ```grep 0  *.tmap | perl -pe 's/Combined\.(\d+)\.gtf\.tmap:/$1\t/' > All.gtf.tmap```
- Get read count data over first exon of TCF4 transcript ENST00000544241.6
    - ```SELECT FPKM from Cufflinks WHERE transcript_id = 'ENST00000544241.6'```
        - no expression of transcript
    - Use mpileup to get read counts over exon (chr18:55403601-55403997)
    - ENSG00000196628
    - In the GTEx junction data, the intron after this exon is 18_53070749_53070852
        - Running mpileup on everything creates a VERY large file, so I've modified mpileup.sh to run over a specified region for all samples
- Analyse WASP results
    - ```samtools mpileup -d 8000 -f /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa -r chr2:184936178-184936178 /c8000xd3/rnaseq-heath/Mappings/15240/BAM/15240.chr.keep.merged.sorted.bam | ~/LabNotes/Python/CountBases.py 15240```
        - this produces no output because all of the reads in this region have been discarded
        - when I run it on the WASP output using only the non-reference bases, I get the same results as I did for the non-remapped BAM file (50 non-ref, zero ref)
        - this makes sense, since this strain is homozygous for the non-ref base (```bcftools view -r chr2:184936178 /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr2.GRCh38.sort.vcf.gz | cut -f 10```)
        - this also means that the focal SNP is not enough to affect the mapping
        - there are 4 bp indels 48 bp upstream and 61 bp downstream from the focal SNP, which explains why all of these reads are lost
        -these indels are not present an ANY of our samples, so it definitely doesn't make sense to use the entire VCF
        - I need to at least filter out any variants absent from all samples
                    

# Cell-type deconvolution
- there is a nice description in the supplement of the common mind paper on how they did this
    - they used [this](http://web.stanford.edu/group/barres_lab/barreslab_rnaseq.xlsx) 7 cell type mouse dataset from {Zhang:2014bt}
    - I need to convert the MGI Symbols to HUGO symbols
        - The 'Complete List of Human and Mouse Homologs (plain text)' link [here](http://www.informatics.jax.org/homology.shtml) has mouse REFseq Gene IDs and symbols for their homologs
        - Now I just need to figure out a way to get REFSEQ IDs from MGI symbols
        - I can do a batch query [here](http://www.informatics.jax.org/batch) that I think will give me what I need
        - I created MGI table to store search results and MouseHumanHomo to store homology INFO
        - I wrote MouseHomo.py to convert symbols in Reference/barreslab_rnaseq.txt to human symbols (MouseHomo.txt)
            - Of 22458 genes with expression info, this was able to assign human homologs to 16094
            - 542 had screwed up search results (no result or multiple). It's possible that if I used Ensembl or something like that to search I might find homologs for some of these
                - The list of homologs includes 85 previous HUGO symbols, 468  Synonyms, and 54 unmatched Symbols
                - If necessary, I can recover all but the unmatched ones using the [HUGO symbol checker](http://www.genenames.org/cgi-bin/symbol_checker)

# Differential Expression Analyses
## HTSeq-count
- Install:
    - ```pip install HTSeq```
    - run htseq-count.sh
    - wrote SubmitHtseq-count.sh to run on all samples sequentially      

## DESeq2/EdgeR
- Used [SARTools](https://github.com/PF2-pasteur-fr/SARTools) script to run analysis
- Ran on output of HTseq-count (using default settings) on the .chr bam file (excludes rDNA)
- Applied an arbitrary grouping to the analysis so there shouldn't be any DE Genes
- 17198 looks really strange in the DESeq analysis because it has very few reads overlapping coding features. I am rerunning with this sample excluded (17025 is also excluded because it has very few reads mapping)
- report includes normalised counts for each gene for each sample (in Counts/tables/num2vsnum1.complete.txt) that can be used to see if a given SNP affects expression of a gene of interest
- Background table includes all genes expressed above the threshold for inclusion in testing. For EdgeR, this can be easily done by filtering out genes were the DE analysis is set to NA, but for DESeq, I manually set genes with Cooks distances > 0.75 to NA, regardless of expression level. This means that some genes that are expressed are also set to NA. I filtered the list of genes myself by extracting the filtering threshold and filtering out genes below the threshold. This doesn't give exactly the same numbers as the number reported as being removed, but it is in the ballpark. Differences are probably due to rounding errors or something
- This manual setting of p-values to NA was necessary because Cooks distance filtering is not done when the formula includes numeric values, which RIN is (also PCW in the full analysis). This is something to do with comparing numbers between groups, which you can't do with numeric variables. I plotted the counts for a lot of genes with elevated Cooks distances, ans 0.75 seemed to be a reasonable cutoff. It is a bit of a hack though.
- This is now set up to run on the server, analysing all data and data divided by week, using both DESeq and EdgeR.
    - I have been using the anaconda version of R because packages that I installed uing the system version didn't always work when submitted jobs to the cluster. [This](http://maven.psycm.cf.ac.uk/wiki/index.php/Scripting_R_library_installation) should solve these problems. 
    - Generally, nominal p-values are lower with DESeq2 than they are with EdgeR
    - EdgeR misses a bunch of genes that are DE according to DESeq because the coverage cutoff (1 CPM) is much more stringent than the DESeq independent filtering threshold
    - DESeq misses a bunch of genes that are DE according to EdgeR because they have samples with high Cooks distances. However, there's only one gene from the full dataset that this applies to (ENSG00000237973). Visual inspection of the data indicates that this gene is expressed at a low level in most samples but at a high level in a subset, only on of which is female.
- It also does DESeq analyses excluding samples with anomalous sequencing depth
- It also analyses changes over development time, in both sexes, as well as an interaction between them.
    - Initially, the coefficients that this was giving me were nonsense. This was because by default, the results() function of DESeq gives coefficients for dropping the last term from the model (in this case, ReadLength). Why it doesn't default to giving coefficients for dropping the term that is dropped in the reduced model is beyond me.
    - After supplying the name='PCW' argument to results(), I get usable coefficients
    - baseMean is the estimated expression for a sample with mean age and log2foldchange is the change per week.
    - The formula to fit a curve is ```baseMean*2^(log2FoldDiff*(PCW-mean_age))```
    - There should be a way to feed this formula to ```geom_smooth()```, but I just calculated it for each week, and used ```geom_line()```
    - This is now implemented in the Shiny app
- It also runs with sex chromosomes features excluded
    - This produces identical nominal p-values, but much higher FDR q-values.
    - This is because the q values are p-values scaled by the number of tests *divided by* the rank of the test ```(q =p*m/i)```. Having a bunch of extra tests at the top of the lists means that marginal tests are ranked lower and more likely to be deemed significant.
    - If you add extra significant tests, m increases but m/i decreases.
        - when sex chromosomes are included, the threshold for independent filtering (which is set automatically to maximum the number of significant tests at q=0.1) is set to 0.35, resulting in 32788 features being excluded, and 30 DE autosomal Genes
            - for ENSG00000189134.3, ```q=0.0000651*25351/16=0.114```
        - when sex chromosomes are excluded, the threshold is 1.7, resulting in 36436 features being excluded (a larger number despite the fact that the total number of features at the beginning is smaller).
        - for ENSG00000189134.3, ```q=0.0000722*31369/76=0.030```
        - HOWEVER, this smaller number of tests results in lower adjusted FDR values
- It is also supposed to systematically try excluding samples from the PCW14 analysis, but this isn't running
- Code PCA plot by sequencing batch
    - The info about sex is really important for interpreting these plots, so I figured out how to use different shapes for other factors
    - This isn't trivial because the prcomp object that pca creates doesn't have simple eigenvectors
    - ggbiplot takes care of converting it to an object that can be plotted with ggplot, but the release version doesn't make it easy to customise the aesthetic mappings.
    - Need this version of ggbiplot: devtools::install_github("richardjtelford/ggbiplot", ref = "experimental") which simply converts the prcomp object into something that can be plotted with ggplot
- The output of this analysis is enhanced with info about which chromosomes the genes are on, wether the chromosome is autosomal, X or Y, the HUGO symbol, and the gene type
    - I've also added 'psuedoautosomal' as an additional chrType, using [this](http://www.genenames.org/genefamilies/PAR) list of PAR genes
    - I've plotted sex bias across the X chromosome, and the start of the chromosome which contains the larger of the PARs includes a cluster of male-biased genes. Female-biased genes are spread across the chromosome, which several genes with the largest bias clustering around XIST

## biobroom
- This can be used to extract model results from DESeq/EdgeR (I don't think it works for Sleuth)
- ```res <- tidy(DESeq(DESeqDataSetFromMatrix(...)))```
- fitted values (for factors): mutate(res, fitted=baseMean*2^estimate)
    - this is different from the values that I have been reporting from SARtools, which are simply group means
- log fold changes are a bit complicated to extract because one has to sum estimates from two rows, but in the binary case is approximately: ```filter(res, term == 'SexMale') %>% mutate(LFC = estimate*2)```
    - the correct estimate for any contrast is ```filter(res, term == 'SexMale' | term == 'SexFemale') %>% select(gene, term, estimate) %>% spread(term, estimate) %>% mutate(LFC = SexMale-SexFemale)```
    - This is equivalent to ```results(DESeq(DESeqDataSetFromMatrix(...)), contrast=c('Sex', 'Male', 'Female'))
    - The former is a lot more typing, but it's WAY faster
- lfcSE is calculated the same way, while Wald stat, p-value and q-value all duplicated in both rows in the binary
- Where there are three levels, all stats are the equivalent of dropping the level in that row; ie the following are equivalent:
    - ```filter(res, term =='"ReadLength2.x.125.bp"```  
    - ```results(DESeq(DESeqDataSetFromMatrix(...)), contrast=c("ReadLength", "2.x.100.bp", "2.x.75.bp"))```
    - (though I had to run the latter as ```results(DESeq(DESeqDataSetFromMatrix(...)), contrast=c(0,0,0,0,1,0,1))```
- No idea what happens with >3 levels    
- incidentally, raw counts can be extracted as a matrix with ```counts((DESeq(DESeqDataSetFromMatrix(...)))``` and normalised counts with ``counts(DESeq(DESeqDataSetFromMatrix(...)), normalized=TRUE)```

## DEXSeq
- This is a pretty heavy duty analysis, so I'm going to try to set it up on rocks
    - There's some problem with the system version of R interacting with my installation of anaconda
    - I'm just going to try installing R with anaconda and using that:
        - ```conda install -c r r```
    - I wasn't able to install DEXSeq because RcppArmadillo would not compile. I installed it using coda:
        - ```conda install -c rgrout r-rcpparmadillo```
        - When I tried again to install DEXSeq, it asked me if I want to update RcppArmadillo. I selected 'no'
    - I ran this on 6 samples and the cpu time was 8.3 hrs, which is slower than JunctionSeq, HOWEVER, it only required 4.9 GB of memory, making it possible to use 8 cores with a wallclock time of 1.7 hrs and total maxvmem=39.038G (RunDEXSeq.sh.o25298)
        - note that this is also using a fuller model with RIN and ReadLength
    - On 10 samples and the cpu time was 10.5 hrs, and maxvmem was 5.35. On 8 cores wallclock time was 1.9 hrs and total maxvmem=42.8 (RunDEXSeq.sh.o25309)
        - this appears to be scaling exactly the same as JunctionSeq, using slightly more cpu and slightly less memory
    - On 16 samples the cpu time was 28.9 hrs with a maxvmem=6.6. On 8 cores with 8G of memory, wallclock time was 4.7 hrs and total maxvmem=52.820 (RunDEXSeq.sh.o25325).
        
## JunctionSeq
- This is a soupped up version of DEXSeq that also uses info about splice junctions.
- I ran it successfully on Carolina's data last year, and it looked pretty good, though there were problems with rare exons being called as differentially spliced because of reads from preRNA
- I ran QoRT on all of my data to get counts across all features.
    - I had to run this on the original sorted BAM file from Tophat because The filtered versions were causing failures that reported orphaned reads.
    - This reports a failure to determine strand for a few of our samples.
        - between 19% and 57% of reads are classified as ```frFirstStrand```
            - this is opposite to RSeQC, which reports the majority of reads on the second strand
            - there also appears to be inconsistency between QoRTs (-stranded) and htseq-count (reverse)/tophat (fr-secondstrand), though the QoRTs manual says that -stranded=fr-firststrand in Tophat2
        - between 0.6% and 4% of reads are classified as ```frSecondStrand```
        - between 40% and 79% are classified as ```ambig_noGenes``` (presumably mostly intronic)
        - strand failure appears to be based on a  frFirstStrand/frSecondStrand ratio cutoff of ca. 9 (8.4 < cutoff > 9.2 )
            - our samples have ratios from 6.2 to 90, with 5 samples below the cutoff
    - There is some additional issue with 16483 which takes FOREVER (and a ton of memory) to run
    - This outputs counts for genes (labeled A001), exons (labeled with E and consecutive numbers), known junctions (labeled with J and consecutive numbers and novel junctions (labeled with N and consecutive numbers).
    - Novel junctions are numbered up to 431, which must mean that 431 different novel splice sites were found in a single gene. There aren't 421 novel splice sites in the counts file though, so that must be the total across all samples, with far fewer in the consensus.
    - There are a total of 229 novel junctions in 181 different genes, with up to 6 novel junctions per gene.
    - There is also a file from merging the novel junctions called orphanSplices which appears to include 197 junctions that aren't in genes.
- This takes a LONG time to run on the server and requires a LOT of memory.
    - I currently have one run going (called JunctionSeq) that is meant to be a single process, but it says it's using 10. Not sure what's up with that.
    - A second process (called JunctionSeq16) started out running on 16 processors, but I had to drop it down to 6 to get enough memory (30 GB per core; 20 failed)
    - Since these are taking days, I'm trying to run it on small subsets of my data to get a handle on how it scales (single process, #$ -l h_vmem=12G)
        - I'm calling these JuncitonSeq6
        - It ran on 6 samples in 6.5 hrs with maxvmem=5.615G (RunJunctionSeq6.sh.o25182)
        - It ran on 10 samples in 8.7 hrs with maxvmem=6.086G (job=RunJunctionSeq6.sh.o25192)
        - It ran on 16 samples in 19.6 hrs with maxvmem=9.438G (RunJunctionSeq6.sh.o25193)
        - Doing a quick and dirty extrapolation from a quadratic fit to these data, I estimate that the 130 sample job running on 6 cores will finish after 300 years!
        - I'm keeping this job running for now, but I killed the mystery job to free up resources
        - After filtering out features with mean counts < 100, I ran this on 121k features. The job died after 2.9 hrs (=18.6 hrs since it was run on 8 cores). maxvmem was 149G= 18.6 per core. This is somehow not any faster than analysing the full dataset and seems to be using WAY more memory.
        - After filtering out *transcripts* with low counts (TPM<1 in at least 56 samples according to Kallisto), it ran on 15 samples using 6 cores in 1.7 hrs (CPU time = 5.3 hrs). maxvmem was 128.495G=22G per core (RunJunctionSeqFilter.sh.o26263)
           - this is a significant improvement (4 fold) in CPU over analysing all transcripts, but using twice as much memory for some reason
               - I think this was accidentally run on 8 cores (need to modify this in the R script as well as the bash script)
               - this means that maxvmem was 128.495G/8=16G per core
               - depending on how this scales with sample size, it may still take 75 years to run on the complete dataset
           - On 32 samples, this ran in 8 hrs on 4 cores (CPU time = 22 hrs), maxvmem was 82.567G = 20G per core
           - this is actually slightly less memory than the 15 sample case, but 4x the CPU time. If this scaling holds, I will need 22 * 4 * 4 = 352 CPU hours to run the complete dataset. That's 2 days on 8 cores with 24G per core. Should be doable!

## Kallisto
- [This](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0862-3) paper recommends filtering out *transcripts* with low counts, not counting bins
    - This approach is endorsed in the [DEXSeq manual](https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.pdf) 
- [Kallisto](https://pachterlab.github.io/kallisto/) can be used to estimate transcript abundance. The output from it can also be used as input for [Sleuth](https://pachterlab.github.io/sleuth/about) to test for differential transcript abundance
    - this is not the same as differential transcript usage because it looks at each transcript independently, not relative to the overall expression of the gene
    - If I can't get DEXSeq to run in reasonable time, I will have to settle for differential transcript abundance, as well as looking at junctions in the same way.
- I was able to get Kallisto to run on one of my samples with 100 bootstraps with maxvmem=4.174G and wallclock=3.2 hrs
    - It is possible to run the bootstraps across multiple processors, but there's not much point since I have to run it on each of my 130 samples.
    - I need to borrow the code from MappingPipeline.sh to get the sample names to line up properly with the fastq files
    - All fastq files per sample are combined at this stage, so the code needs to be a bit different.
    - When I run 12107 with --rf-stranded, 33M of 85M reads pseudoaligned. When I run with --fr-stranded, 3M reads are pseudoaligned.
- Kallisto produces a tab-separated output with gene lengths, "effective lengths", estimated counts, and TMP (transcripts per million).
    - I think it makes sense to use TPM to filter out transcripts with low expression.
    - To that end, I'm going to write an R script that produces a data frame with TPM across samples. I can then use this to select a list of features to remove from the GTF file, which can then be used to generate a new GFF file for counting features across exons/junctions (this is for a subset of 117 samples).
    - Of 200k transcripts:
        - 162k with average TPM > 0.01
        - 119K with average TPM > 0.1
        - 50k with average TPM > 1
        - 4.5k with average TPM > 10
        - 136k with TPM > 0.01 in at least 46 samples
        - 111k with TPM > 0.1 in at least 46 samples
        - 50k with TPM > 1 in at least 46 samples
        - 4.7k with TPM > 10 in at least 46 samples
            - for EdgeR, we used CPM > 1 in at least 46 samples
- Using default parameter, Sleuth requires a seemingly infinite amount of RAM
    - I thought this was because of all the bootstrap replicates, so I redid Kallisto wiht 10 reps
    - This turns out to be unnecessary because there's an option to limit the number of reps used by Sleuth HOWEVER, the problem turned out to be that by default, Sleuth uses ALL the processors. Once I limited it to a single core, it ran in reasonable time with reasonable memory   
## Differential splicing
- The most direct way to evaluate differences in splicing is to use the spliced reads from the mapping.
- These can be used for true differential splicing (changes in spliced reads *relative* to reads mapping to the gene) using DEXseq or they can simply be tested for differences in abundance, using eg; DESeq
- The latter approach is a bit statistically dodgy, because junctions within a gene are not independent of each other (ie; if a gene is DE, then all of the junctions that make it up should be DE). However, it has the advantage that it hopefully won't take 300 years to run. Any DE junctions that are not part of a DE gene can be investigated on a case-by-case basis, or I can come up with some systematic approach to dealing with them
- Tophat produces a bed file with splice junctions and the number of reads that support them (229,546 junctions)
- There is also a file with counts over known splice sites output by QoRT, as well as a file that includes novel splices, which could be extracted with ```grep ':N'``` (377,531 + 229 novel)
    - This includes a lot of zero counts that aren't in Tophat output, but they are all matched up between samples and assigned to genes, which makes life a bit easier. The zero counts are just going to get screened out by DESeq anyway.


## Cufflinks
- Ran Cuffmerge.sh to get combined gtf, followed by Cuffquant.sh to get FPKM values
- Cuffquant produces a binary .cxb file. I'll need to run this on everything, then run cuffnorm to get FPKM values

# SNP calling
- This is needed to ensure that the samples used for RNAseq match the ones used for genotyping
- In theory, SNPs called from the RNAseq mappings should match the SNPs in the imputed VCF
- I don't need a huge number of SNPs for this analysis, so I extracted the mappings to chr22:
    - ```qsub ~/LabNotes/SubmissionScripts/DivideBAM.sh 22```
- Call SNPs using samtools mpileup and bcftools call
    - ```find /c8000xd3/rnaseq-heath/Mappings/ -name *chr22.bam | xargs qsub ~/LabNotes/SubmissionScripts/CallSNPs.sh```
- Change names so that they match the ones in the reference VCF
    - ```find /c8000xd3/rnaseq-heath/Mappings/ -name *chr22.bam | perl -pe 's/\/c8000xd3\/rnaseq-heath\/Mappings\/([\w-]+)\/BAM\/Chromosomes\/[\w-]+.chr22.bam/$1/' 
> /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/SampleNames.txt```
    - ```bcftools reheader -s /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/SampleNames.txt -o /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/chr22.renamed.bcf /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/chr22.raw.bcf```
    - ```bcftools index /c8000xd3/rnaseq-heath/Genotypes/SNPcalls/chr22.renamed.bcf```
- Use bcftools gtcheck to confirm that samples match:
    - ```bash ~/LabNotes/SubmissionScripts/SubmitGTcheck.sh```

# [ASEReadCounter](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php)
- This tool
 from GATK seems like it's probably the best bet to extract AS read counts (tho it looks like WASP has a tool now also?)
- Unfortunately, GATK it awful picky about input BAM files. I need to run [ValidateSAM](https://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile) from Picard and fix all the issues with them before proceeding
    - see info [here](https://software.broadinstitute.org/gatk/documentation/article?id=7571)
    - MISSING_READ_GROUP/RECORD_MISSING_READ_GROUP
        - this info can be added with [AddOrReplaceReadGroups](https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups) from Picard
        - according to [this](http://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups) extremely confusing explanation, RGID and RGPU are somehow different, but each is unique for each sample/lane
        - I have opted to use the flow cell ID, lane, and barcode to fill both of these fields
            - ```find /c8000xd3/databank/foetal-rna/ -name $(grep 15240 ~/LabNotes/sequences.txt | head -1 | cut -f 1)* | xargs zcat | head -1 | cut -d: -f 3,4,10 | perl -pe 's/:/./g'```
        - RGLB is 'lib1' (or lib2 in the case of 15533-2) and RGPL is illumina
        - RGSM is the sample name
    - MATE_NOT_FOUND
        - this is triggered if the read is flagged as paired in sequencing ([flag 1](https://broadinstitute.github.io/picard/explain-flags.html)), but mate not present in the bam file
        - this is the case for me whenever the mate is not mapped ([flag 8](https://broadinstitute.github.io/picard/explain-flags.html))
        - [FixMateInformation](https://broadinstitute.github.io/picard/command-line-overview.html#FixMateInformation) from Picard can supposedly fix this (I think. the documentation just says 'fix if needed'), but it doesn't appear to be doing anything to my mappings
        - I'm actually really confused by this because the output file is 1 GB bigger than the input, but ValidateSAM gives the same number of missing mates. I guess it's done SOMETHING to the non-missing ones?
            - Tophat produces a second bam file called unmapped.bam. I'm hoping that if I combine this with accempted_hits.bam, all will be right with the world
            - It looks like [MergeBamAlignment](https://broadinstitute.github.io/picard/command-line-overview.html#MergeBamAlignment) from picard is designed to do just this
            - This is a nightmare. There's all manner of problems with the unmapped.bam file form tophat that cause picard to choke
            - fortunately, I'm not the first to encounter this problem, and I've found a tool called [tophat-recondition](https://github.com/cbrueffer/tophat-recondition) that is meant to solve these problems
            - [AddOrReplaceReadGroups](https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups) is now complaining about the sort order. I am going to try [SortSam](https://broadinstitute.github.io/picard/command-line-overview.html#SortSam) to sort the unmpped reads by queryname and see if that helps
            - This seems to work, but I should also run it on the accepted_hits because [MergeBamAlignment](https://broadinstitute.github.io/picard/command-line-overview.html#MergeBamAlignment) complains about it
            - Alternatively, I'm going to try coordinate-sorting the unmapped file (which seems like a weird thing to do) to see if that works
                - coordinate sorting the unmapped file before adding read groups does work, but it doesn't make the warning (Exception merging bam alignment - attempting to sort aligned reads and try again: Underlying iterator is not queryname sorted) doesn't go away
            - [MergeBamAlignment](https://broadinstitute.github.io/picard/command-line-overview.html#MergeBamAlignment) STILL doesn't work. Now it's complaining about identical Program Record ID's in the headers from the mapped and unmapped files
            - This is because both have ID:Tophat in the @PG line of the header. I think I have managed to solve this by adding a line to [tophat-recondition](https://github.com/cbrueffer/tophat-recondition) to modify the PGID of the unmapped file to be "Tophat-unmapped". This is running right now, so we'll see how it goes
            - I also had to modify [tophat-recondition](https://github.com/cbrueffer/tophat-recondition) to set the mate_reverse_strand flag
            - ValidateSAM is still complaining about the PNEXT for unmapped reads being set to one, rather than the position of the mapped read and about RNEXT/PNEXT of the mapped read being set to */0, rather than the coordinates of the mapped read
                - this is because tophat-recondition set the RNAME/POS of unmapped reads to RNAME/POS of the mapped mate, but does not change RNEXT/PNEXT of the mapped mate
                - changeing PNEXT of the unmapped reads is an easy fix by tweaking tophat-recondition, but that script doesn't modify accepted_hits, so fixing the mapped mate is not trivial
                - I'm trying MergeBAM from picard on the latest output of tophat-recondition to see if it can cope with this, but I'm not optomistic
                    - when reads are queryname sorted, MergeBAM gives an error: "picard.PicardException: Second read from pair not found in unmapped bam: HISEQ:258:C6UARANXX:4:1101:10000:28046, HISEQ:258:C6UARANXX:4:1101:10000:64326"
                    - when reads are coordinate sorted, MergeBAM first does a queryname sort, then gives error
                - It's also possible that ASEReadCounter will ignore these problems that ValidateSAM is flagging
                - I've modified tophat-recondition.py to make the necessary changes to accepted_hits.bam
                    - if mapped.mate_unmapped == true, mapped.pnext = mapped.pos, mapped.rnext = mapped.tid
             - No error reported!
             - I've incorporated all the steps to do this into a script called Tophat2GATK.sh       
        - WASP somehow fixes this missing mate problem, but creates a new problem where reads with properly mapped mates have the mates removed
        - I suspect that this is due to poor mate read quality
    
        	| BAM_file	                                 | MATE_NOT_FOUND |
        	| ------------------------------------------ | -------------- |
        	| accepted_hits.bam	                         | 4357912        |
        	| 15240.chr.bam		                         | 4319245        |
        	| 15240.chr.sort_RG_mateFix.bam		         | 4319245        |
        	| 15240.chr.nonref.merged.bam		         | 121069         |
        	| 15240.chr.nonref.merged.sorted.bam		 | 121069         |
        	| 15240.chr.nonref.merged.sorted_mateFix.bam | 121069         |
        	| 15240.chr.nonref.merged.sorted_RG.bam      | 121069         |

    - SAMException: Value was put into PairInfoMap more than once
        - somewhat ironically, the WASP deduplication script duplicated the entry for one of the reads
        - hopefully this is fixed with the new version of WASP
    - I'm debating the proper order of doing things here. I managed to get the output from Tophat to pass VerifyBAM, but that was without WASP remapping. I can run WASP on this GATK compatible file, but it requires tophat so I think I'm just going to reintroduce all the same problems. Probably the thing to do is run WASP as the first step after tophat, but what about the remove duplicates step?
        - Sorting and adding read groups appears to be enough to get the mapping to pass ValidateSAM
        - I've written an R script to combine the counts from individual mappings. This doesn't give exactly the same counts as passing multiple fastq files to Tophat, but it's close (a total of 159 features with discrepant counts, with a total of 206 differences). I think I'm not going to worry about it.
    - This seems to work fine BUT, it required filtering out duplicated sites from the vcf files (ie; sites with >2 genotypes)
        - I've updated ProcessVCF.sh to filter out these sites, but the MD5 files used for WASP were made from vcf files that include these sites. This is probably preferable, but it needs to be noted
