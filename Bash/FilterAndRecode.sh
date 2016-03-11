~/bin/plink --bfile FB_Merged --chr $1 --recode vcf --out FB_Merged_chr$1
~/BTSync/Code/Python/checkVCFmod.py -r ~/Documents/src/checkVCF-20140116/hs37d5.fa -o FB_Merged_chr$1 FB_Merged_chr$1.vcf
egrep 'Ambiguous|Ref not found' FB_Merged_chr$1.check.ref |cut -f 3 >FB_Merged_chr$1_exclude.txt
cut -f 3 FB_Merged_chr$1.check.dup >>FB_Merged_chr$1_exclude.txt
grep 'Strand switch' FB_Merged_chr$1.check.ref | cut -f 3 > FB_Merged_chr$1_flip.txt
plink --bfile FB_Merged --chr $1 --exclude FB_Merged_chr$1_exclude.txt --flip FB_Merged_chr$1_flip.txt --recode vcf --out FB_Merged_chr$1_flip
vcf-sort FB_Merged_chr$1_flip.vcf |bgzip -c > Filtered/FB_Merged_chr$1_flip.vcf.gz
~/BTSync/Code/Python/checkVCFmod.py -r ~/Documents/src/checkVCF-20140116/hs37d5.fa -o FB_Merged_chr$1_flip FB_Merged_chr$1_flip.vcf
