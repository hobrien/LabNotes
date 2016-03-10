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
    - ```find /c8000xd3/databank/foetal-rna/Exeter\ sequencing/ -name /*fastqc/* -print0 | xargs -n 1 bash SubmissionScripts/cp.sh```
    - FastQC.sh:
        - ```~/src/FastQC/fastqc --outdir=/home/heath/FastQC "$@"```
    - cp.sh:
        - ```cp "$@" >FastQC/```