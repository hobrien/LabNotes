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
            - ```grep MismatchRefBase FB_Merged_chr1.report.check.ref | cut -d : -f 2 | sort -b > FB_Merged_chr1_flip_pos.txt```
        - these need to be compared to the VCF to identify the SNPs at those positions
            - ```grep -v '#' FB_Merged_chr1.vcf | cut -f 2,3 | sort -b |join - FB_Merged_chr1_flip_pos.txt | cut -f 2 -d ' ' | sort -n >FB_Merged_chr1_flip.txt```
            - this didn't work because there are two SNPs at position 11854476
                - VG01S2022 and rs1801131
        - flip SNPs
            - ```plink --bfile FB_Merged --chr 1 --flip FB_Merged_chr1_flip.txt --recode vcf --out FB_Merged_chr1_flip```
        - rerun checkVCF.py
            - ```checkVCF.py -r ~/Documents/src/checkVCF-20140116/hs37d5.fa -o FB_Merged_chr1_flip FB_Merged_chr1_flip.vcf ```    