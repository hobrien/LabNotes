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
  * [Transcript Identification](#transcript-identification)
    * [Cufflinks](#cufflinks)
    * [Cuffcompare](#cuffcompare)
    * [Make DB of GencodeGTF:](#make-db-of-gencodegtf)
  * [Expression analysis](#expression-analysis)
  * [Cell\-type deconvolution](#cell-type-deconvolution)

Created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc.go)


#Create VCF with imputed SNPs
##Prepare files for imputation
###PLINK
- Followed the instructions [here](https://imputationserver.sph.umich.edu/start.html#!pages/help) to convert binary PLINK (.bed) files to input for the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/start.html)
    - I checked some of the SNPs on the [USCS Genome Browser](http://genome-euro.ucsc.edu/index.html) and it looks like it uses GRCh37/hg19, but the SNPs appear to be on the opposite strand
    - Nick is going to check and try to figure out what is up with the chr0 SNPs (1181 of them)
    - need to separate sequences by chromosome first
        - ```echo {1..25} |xargs -n 1 -P 8 -I % ~/bin/plink --bfile FB_Merged --chr % --recode vcf --out FB_Merged_chr%```
        - this gave warnings about hets on the haploid chromosomes and non-missing data on the Y.
        - The 103 males all have between 26 and 96 heterozygous SNPs on chr23 and between 3 and 17 on the Y
        
    - sort and compress 
        - ```echo {1..25} |xargs -n 1 -P 8 bash SortAndCompress.sh```
        - SortAndCompress.sh: ```vcf-sort FB_Merged_chr$1.vcf |bgzip -c > FB_Merged_chr$1.vcf.gz```
        - I did this because I couldn't figure out how to use xargs with a pipe

###checkVCF
- Run checkVCF.py
    - ```checkVCF.py -r ~/Documents/src/checkVCF-20140116/hs37d5.fa -o FB_Merged_chr1 FB_Merged_chr1.vcf ```
    - this gave a lot of wordy warnings which the server help says to ignore, but the VCFs appear to be valid
    - I didn't bother to run it on all the files
    
- I was able to upload the files for the autosomes, but the sex chromosomes are giving me trouble
    - Nick says not to worry about them

- Turns out the the flipped SNPs are big problem for imputation. They need to be fixed.
    - checkVCF.py produces a file called XXX.report.check.ref with a list of flipped SNP positions
    - there are also two SNPs at position 11854476
        - VG01S2022 and rs1801131
    - I modified the script (checkVCFmod.py) to report the SNP id and the kind of switch
        - options are:
            - 'Allele switch'
            - 'Strand switch'
            - 'Strand switch and allele switch' (requires both)
            - 'Ambiguous' (either switch would produce correct ref)
            - 'Ref not found' (neither the alt allele or the rev of either matches the ref)
        - Chr1 has 7741 allele switches, 21446 strand switches, 7905 with both, 127 ambiguous and 19 where the ref is not found
            - NB: grep 'Strand switch' gets all strand switches, including those with both
    - Flip strand:
        - ```grep 'Strand switch' FB_Merged_chr1_mod.check.ref | cut -f 3 > FB_Merged_chr1_flip.txt```
        - ```plink --bfile FB_Merged --chr 1 --flip FB_Merged_chr1_flip.txt --recode vcf --out FB_Merged_chr1_flip```          

    - exclude ambiguous SNPs, SNPs missing reference and duplicates
        - I modified checkVCF.py to include the name of the duplicate SNPs, so this will make a list
            - ```egrep 'Ambiguous|Ref not found' FB_Merged_chr1_mod.check.ref |cut -f 3 >FB_Merged_chr1_exclude.txt```
            - ``` cut -f 3 FB_Merged_chr1_mod.check.dup >>FB_Merged_chr1_exclude.txt```
        - This will make a VCF without these SNPs
            - ```plink --bfile FB_Merged --chr 1 --exclude FB_Merged_chr1_exclude.txt --flip FB_Merged_chr1_flip.txt --recode vcf --out FB_Merged_chr1_flip_filter```
    - Fixing the allele flips doesn't appear to be so straightford, but they just say INFO in the imputation server output, not FILTER, so hopefully they are ok
    - repeat sorting and compressing on flipped/filtered vcf
        - ```vcf-sort FB_Merged_chr1_flip_filter.vcf |bgzip -c > FB_Merged_chr1_filtered.vcf.gz```
    - I'm not yet sure if this will work or not, but I might as well figure out a way to apply this to all chromosomes
        - ```echo {1..22} |xargs -n 1 -P 4 bash ~/BTSync/FetalRNAseq/LabNotes/Bash/FilterAndRecode.sh```
- Run imputation (see ScreenShots/Imputation2.png)
    - Imputed data are in BTSync/FetalRNAseq/ImputedGenotypes/Raw_output/    

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
        
        - ``` find Uncompressed/ -name summary.txt | xargs perl -pe 's/\..*//' >>summary.txt```
    - Extract the number of reads for each SampleID
        - ```find Uncompressed/ -name fastqc_data.txt | xargs grep 'Total Sequences' | grep -v trimmed | perl -pe 's/(\.sanfastq)?_fastqc\/fastqc_data\.txt\:Total Sequences//' | perl -pe 's/\.\///' > seq_lengths.txt```
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
            - That would be because I selected the unphased output option
                            
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

##Cufflinks
- Ran Cuffmerge.sh to get combined gtf, followed by Cuffquant.sh to get FPKM values
- Cuffquant produces a binary .cxb file. I'll need to run this on everything, then run cuffnorm to get FPKM values
