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
    - Concatenate and filter imputed data    
        - Imputed data are in BTSync/FetalRNAseq/ImputedGenotypes/Raw_output/
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
            - ```bcftools reheader -h temp.head -o subset2.vcf.gz subset.vcf.gz```
            - ```rm temp.head```
    - replace genotyping IDs with Brain Bank IDs
        - ```python ~/BTsync/FetalRNAseq/LabNotes/Python/ChangeSampleID.py All_chromosomes.vcf.gz```
        - ```bcftools reheader -h temp.head -o All_chromosomes_relabelled.vcf.gz All_chromosomes.vcf.gz```
        - ```rm temp.head```

    - Add SNP IDs to imputed VCFs
        - SNP DB [Schema](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp144.sql) and [data](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp144.txt.gz) downloaded from the [USCSC Genome Browser](https://genome.ucsc.edu/) and imported into FetalRNAseq mySQL DB
        - chromosome (field 1), position (2) and ref (4) and alt (5) alleles can be used to uniquely identify each SNP (I hope!)
            - there can be several SNPs at the same position
        - these correspond to chrom, chromStart, refNCBI and observed in snp144
                                
        - ```grep -v '^#' FB_Merged_chr22.vcf | cut -f 1,2,3 |python ../LabNotes/Python/AddSNPID.py | wc -l```
        
        
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

- Trim Adaptors from reads using cutadapt
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

- Start mapping reads with tophat
    - follow the instructions [here](http://www.illumina.com/documents/products/technotes/RNASeqAnalysisTopHat.pdf) to get started with tophat
        - ```wget --ftp-user=igenome --ftp-password=G3nom3s4u ftp://ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz```
        - ```tar -xzf Homo_sapiens_UCSC_hg19.tar.gz```
        - ```tophat2 --GTF /home/heath/Index/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf --library-type firststrand --mate-inner-dist 300 --mate-std-dev 50 --num-threads 8 --output-dir Mappings/15533/ /home/heath/Index/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome /home/heath/Trimmed/15533_TGACCA_L007_R1_001_trimmed.fastq.gz /home/heath/Trimmed/15533_TGACCA_L007_R2_001_trimmed.fastq.gz```
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
    - Try to run RNAseq-specific QC:
        - Eilis used something called RNA-SeQ. She went thru several steps to format the mapping before this would work 
            - ```java -jar ~/Documents/src/picard-tools-1.119/CreateSequenceDictionary.jar R=genome.fa O=genome.bam```
            - ```java -jar ~/Documents/src/picard-tools-1.119/AddOrReplaceReadGroups.jar I=accepted_hits.bam O=accepted_hits2.bam RGLB="totalRNA" RGPL="Illumina" RGPU="1" RGSM="ID_15533"```
            - ``` samtools index accepted_hits2.bam```
            - ```java -jar ~/Documents/src/picard-tools-1.119/ReorderSam.jar I=accepted_hits2.bam O=accepted_hits3.bam R=~/BTSync/FetalRNAseq/Reference/genome.fa```
            - ```samtools index accepted_hits3.bam```
            - ```java -jar ~/Documents/src/RNA-SeQC_v1.1.8.jar -n 1000 -s "ID_15533|accepted_hits3.bam|Test" -t ~/BTSync/FetalRNAseq/Reference/gencode.v19.annotation.gtf -r ~/BTSync/FetalRNAseq/Reference/genome.fa -o RNAseQC -gc ~/BTSync/FetalRNAseq/Reference/gencode.v7.gc.txt```
                - still not working 
        - Run [CollectRnaSeqMetrics]( http://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics) from [Picard Tools](http://broadinstitute.github.io/picard/)
            
            - had to install a more recent version of java on rocks
            
            - downloaded a correctly formatted [file](https://gist.github.com/slowkow/b11c28796508f03cdf4b) with rRNA intervals     
                - I'm pretty sure this stat is not correct. I don't know if it's a problem with this file or with the command I used
            - I had to build a dictionary for the reference Seq
                - ```/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=/home/heath/Ref/hg19.fa O=/home/heath/Ref/hg19.dict```    
            - A (sort of) explanation of the output is [here](https://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics)
        - Run [RSeQC](http://rseqc.sourceforge.net)
            - Downloaded a [BED](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens) to Reference
            - Get Chromosome sizes
                - ```bash Lab_notes/Bash/fetchChromSizes hg19 >Reference/hg19.chrom.sizes```
            - ```bam_stat.py -i ~/Documents/Mappings/15533_300/accepted_hits.bam```
            - ```inner_distance.py -i accepted_hits.bam -o ~/BTSync/FetalRNAseq/BamQC/15533_1000/15533_1000 -r ~/BTSync/FetalRNAseq/Reference/hg19_RefSeq.bed -u 5000 -s 50 > /dev/null```    
                    

- Analyse expressed SNPs
    - Run mpileup on SNPs from grant:
        - ```cat ~/BTSync/FetalRNAseq/Info/ExpressedSNPs.txt | python ~/BTSync/FetalRNAseq/LabNotes/Python/GetSNPpos.py | xargs -n 1 -I % samtools mpileup -d 8000 -f ~/BTSync/FetalRNAseq/Reference/genome.fa -r % -ABQ 0 accepted_hits.bam |python ~/BTSync/FetalRNAseq/LabNotes/Python/CountBases.py ```        
