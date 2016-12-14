  * [Create VCF with imputed SNPs](#create-vcf-with-imputed-snps)
    * [Prepare files for imputation](#prepare-files-for-imputation)
      * [PLINK](#plink)
      * [checkVCF](#checkvcf)
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
    * [Mapping QC](#mapping-qc)
    * [Try to run RNAseq\-specific QC:](#try-to-run-rnaseq-specific-qc)
      * [RNA\-SeQ](#rna-seq)
      * [CollectRnaSeqMetrics](#collectrnaseqmetrics)
      * [RSeQC](#rseqc)
  * [Allele\-specific alignment](#allele-specific-alignment)
    * [I am going to need to do a liftover before running this because the vcf is for hg19](#i-am-going-to-need-to-do-a-liftover-before-running-this-because-the-vcf-is-for-hg19)
      * [LiftoverVcf (picard\-tools)](#liftovervcf-picard-tools)
      * [CrossMap](#crossmap)
    * [WASP](#wasp)
    * [Clip overlapping reads](#clip-overlapping-reads)
  * [Transcript Identification](#transcript-identification)
    * [Cufflinks](#cufflinks)
    * [Cuffcompare](#cuffcompare)
    * [Make DB of GencodeGTF:](#make-db-of-gencodegtf)
    * [Run Cufflinks on SRA data from Liu:2016ji (GSE71315)](#run-cufflinks-on-sra-data-from-liu2016ji-gse71315)
  * [Expression analysis](#expression-analysis)
  * [Cell\-type deconvolution](#cell-type-deconvolution)
    * [HTSeq\-count](#htseq-count)
    * [DESeq2](#deseq2)
      * [DEXSeq](#dexseq)
    * [Cufflinks](#cufflinks-1)
  * [SNP calling](#snp-calling)
  * [ASEReadCounter](#asereadcounter)

Created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc.go)


#Create VCF with imputed SNPs
##Prepare files for imputation
- All steps needed for this analysis are run by the check-bim script on [this](http://www.well.ox.ac.uk/~wrayner/tools/) site (v4.2.5, with a few modifications for additional processing steps)
    - changed output from bed to VCF
    - added sort and compression steps
    - final output written to subfolder, which makes upload to [Michigan imputation server](https://imputationserver.sph.umich.edu/start.html) easier
- Analysis steps:
    - ```cd BTSync/FetalRNAseq/Genome-wide\ genotyping/```
    - ```plink --bfile FB_Merged --freq```
    - ```mkdir Imputation3```
    - ```cd Imputation3```
    - ```perl ../LabNotes/Perl/HRC-1000G-check-bim.pl -b ../FB_Merged.bim -f ../plink.frq -r ../../Ref/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h```
    - ```bash Run-plink.sh```
- Run imputation using HRC r1.1 2016 and Eagle v2.3 for phasing with Mixed population for QC (see ScreenShots/Imputation3.png)
    - files are uploaded from BTSync/FetalRNAseq/Genome-wide\ genotyping/Imputation3
    - Imputed data are in BTSync/FetalRNAseq/ImputedGenotypes/Imputation3 
        
##Concatenate and filter imputed data
###bcftools
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

##Add SNP IDs to imputed VCFs
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
        
#Prepare Sequencing data
##Sequence data QC
###FastQC        
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

##Trimming sequences
###Trim Adaptors from reads using cutadapt
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

###Use trim_galore cutadapt wrapper
- For some obscure reason, the reads are becoming unpaired after trimming. I get very high pair congruence when mapping the raw data with Tophat, but very low congruence with the trimmed data
- I've also decided that it's a good idea to do some quality trimming because unlike Star, Tophat does not do read clipping
- [trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore) is a perl wrapper for cutadapt that makes it a little more flexible and easy to use. It also solves this problem
    - After doing the mapping with quality trimming and adapter trimming, everything looks good except that there's a big spike in SNP frequency at the 5' end of the reads.
    - I am going to solve this by clipping 5 bp from the 5' end of each read in addition to the quality and adapter trimming:
        - ```find /c8000xd3/databank/foetal-rna/Exeter_sequencing/15533_Oct_2014/ -name *.fastq.gz | xargs qsub SubmissionScripts/Trim.sh```

#Mapping            
##Start mapping reads with tophat
- follow the instructions [here](http://www.illumina.com/documents/products/technotes/RNASeqAnalysisTopHat.pdf) to get started with tophat
    - ```wget --ftp-user=igenome --ftp-password=G3nom3s4u ftp://ftp.illumina.com/Homo_sapiens/NCBI/GRCh38Decoy/Homo_sapiens_NCBI_GRCh38Decoy.tar.gz```
    - ```tar -xzf Homo_sapiens_UCSC_hg19.tar.gz```
    - ```tophat --keep-fasta-order --transcriptome-index /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.inx --library-type fr-secondstrand --mate-inner-dist 150  --mate-std-dev 50 --num-threads 8 --output-dir /home/heath/Mappings/15533_150 /home/heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome /home/heath/Trimmed/15533_TGACCA_L007_R1_001_trimmed.fastq.gz /home/heath/Trimmed/15533_TGACCA_L007_R2_001_trimmed.fastq.gz```

##Mapping QC
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

##Try to run RNAseq-specific QC:
### RNA-SeQ
- Eilis used something called RNA-SeQ. She went thru several steps to format the mapping before this would work 
    - ```java -jar ~/Documents/src/picard-tools-1.119/CreateSequenceDictionary.jar R=genome.fa O=genome.bam```
    - ```java -jar ~/Documents/src/picard-tools-1.119/AddOrReplaceReadGroups.jar I=accepted_hits.bam O=accepted_hits2.bam RGLB="totalRNA" RGPL="Illumina" RGPU="1" RGSM="ID_15533"```
    - ``` samtools index accepted_hits2.bam```
    - ```java -jar ~/Documents/src/picard-tools-1.119/ReorderSam.jar I=accepted_hits2.bam O=accepted_hits3.bam R=~/BTSync/FetalRNAseq/Reference/genome.fa```
    - ```samtools index accepted_hits3.bam```
    - ```java -jar ~/Documents/src/RNA-SeQC_v1.1.8.jar -n 1000 -s "ID_15533|accepted_hits3.bam|Test" -t ~/BTSync/FetalRNAseq/Reference/gencode.v19.annotation.gtf -r ~/BTSync/FetalRNAseq/Reference/genome.fa -o RNAseQC -gc ~/BTSync/FetalRNAseq/Reference/gencode.v7.gc.txt```
        - still not working 
        
###CollectRnaSeqMetrics
- Run [CollectRnaSeqMetrics]( http://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics) from [Picard Tools](http://broadinstitute.github.io/picard/)
            
    - had to install a more recent version of java on rocks
            
    - downloaded a correctly formatted [file](https://gist.github.com/slowkow/b11c28796508f03cdf4b) with rRNA intervals     
        - I'm pretty sure this stat is not correct. I don't know if it's a problem with this file or with the command I used
    - I had to build a dictionary for the reference Seq
        - ```/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=/home/heath/Ref/hg19.fa O=/home/heath/Ref/hg19.dict```    
    - A (sort of) explanation of the output is [here](https://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics)

###RSeQC
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

#Allele-specific alignment

##I am going to need to do a liftover before running this because the vcf is for hg19
###LiftoverVcf (picard-tools)    
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

###CrossMap
- fortunately, it looks like something called [CrossMap](http://crossmap.sourceforge.net) is going to come to my rescue. 
    - It works fine on 100,000 lines. We'll see how it does on all of them
    - The only problem is that it produces an uncompressed VCF, so I'll have to run bgzip after it finishes    

##WASP
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
            - I changed the names of the samples in the VCF headers to match the samples
                - ```cut -f 1 ~/LabNotes/VCFindex.txt | cut -d'/' -f 1 | cut -d'-' -f 1 > ~/LabNotes/VCFindex2.txt```
                - ```for chr in {1..22}; do bcftools reheader -s ~/LabNotes/VCFindex2.txt -o /c8000xd3/rnaseq-heath/Genotypes/Imputation2/Sorted/chr$chr.GRCh38.sort.vcf.gz /c8000xd3/rnaseq-heath/Genotypes/Imputation2/chr$chr.GRCh38.sort.vcf.gz; done```
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
                - ```bcftools filter -e 'POS>=89346516' chr22.GRCh38.vcf.gz > chr22.GRCh38.filter.vcf.gz```
                 
##Clip overlapping reads
- the tool of choice for this appears to be [clipOverlap](http://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap) from [bamUtil](http://genome.sph.umich.edu/wiki/BamUtil)
- this clips the read with the lowest quality score, which isn't as good as comparing the bases of the overlapping reads, but is probably good enough
- this turned out to be a huge pain to install, until [bioconda](https://bioconda.github.io/index.html#setup) saved the day
    - ```conda config --add channels r```
    - ```conda config --add channels bioconda```
    - ```conda install bamutil```
  
                            
#Transcript Identification
##Cufflinks
- Run [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks)
    - Cufflinks is running on 15533, but it's taken 3 days so far and no indication of progress
    - Cufflinks finished after 1 day when I ran it on 15468, which isn't a much smaller dataset, so I'm not sure what's going on with that
    - Cufflinks includes an option (-M) to mask reads matching the specified features (ie; rDNA), but it's probably easier to just use the BAM files that already exclude these reads
    - I think it would also be helpful to get rid of all of the non-chromosome sequences.
        - FilterBAM.sh should do this

##Cuffcompare
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

##Make DB of GencodeGTF:
    - ```cat  gencode.v24.chr_patch_hapl_scaff.annotation.gtf |python ~/BTSync/FetalRNAseq/LabNotes/Python/GTF2CSV.py > ~/BTSync/FetalRNAseq/Reference/gencode.v24.chr_patch_hapl_scaff.annotation.csv```
    - ```mv features.csv ~/BTSync/FetalRNAseq/Reference/gencode.v24.chr_patch_hapl_scaff.features.csv```
    - ```mysql db_name < ~/BTSync/FetalRNAseq/LabNotes/SQL/refGTF.sql```
    - ```SELECT DISTINCT feature FROM GencodeFeatures```
        - transcript_support_level, exon_number, level, gene_type, transcript_type, gene_name, transcript_name, tag, ccdsid, exon_id, gene_id, protein_id, transcript_id, gene_status, transcript_status, havana_gene, havana_transcript, ont        

##Run Cufflinks on SRA data from Liu:2016ji (GSE71315)
- downloaded sra toolkit and ran fastq-dump (GetSRA.sh)
- hopefully I get similar results to what I have for our data
- I've downloaded all bulk tissue riboZero RNAseq datasets and Tophat is running on one
- The mapping QC looks VERY similar to our data, but a lot less reads, resulting in 10+ fold fewer novel junctions
    - This suggests that using a higher coverage threshold is the answer to cleaning up the novel cufflinks data
      
#Expression analysis
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
                    

#Cell-type deconvolution
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

##HTSeq-count
- Install:
    - ```pip install HTSeq```
    - run htseq-count.sh
    - wrote SubmitHtseq-count.sh to run on all samples sequentially      

##DESeq2
- Used [SARTools](https://github.com/PF2-pasteur-fr/SARTools) script to run analysis
- Ran on output of HTseq-count (using default settings) on the .chr bam file (excludes rDNA)
- Applied an arbitrary grouping to the analysis so there shouldn't be any DE Genes
- 17198 looks really strange in the DESeq analysis because it has very few reads overlapping coding features. I am rerunning with this sample excluded (17025 is also excluded because it has very few reads mapping)
- report includes normalised counts for each gene for each sample (in Counts/tables/num2vsnum1.complete.txt) that can be used to see if a given SNP affects expression of a gene of interest

###DEXSeq
- This is a pretty heavy duty analysis, so I'm going to try to set it up on rocks
    - There's some problem with the system version of R interacting with my installation of anaconda
    - I'm just going to try installing R with anaconda and using that:
        - ```conda install -c r r```
    - I wasn't able to install DEXSeq because RcppArmadillo would not compile. I installed it using coda:
        - ```conda install -c rgrout r-rcpparmadillo```
        - When I tried again to install DEXSeq, it asked me if I want to update RcppArmadillo. I selected 'no'
          
##Cufflinks
- Ran Cuffmerge.sh to get combined gtf, followed by Cuffquant.sh to get FPKM values
- Cuffquant produces a binary .cxb file. I'll need to run this on everything, then run cuffnorm to get FPKM values

#SNP calling
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

#[ASEReadCounter](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php)
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