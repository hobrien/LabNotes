#!/usr/bin/env python

from __future__ import print_function
import sys
import warnings
from pysam import VariantFile

"""
For SNP calling output, I would like counts of the following genotypes for each SNP:
0/0
0/1 + 1/0
1/1
0/. + ./0
1/. + ./1
./.

- I need to make sure the the delimiter is \ and not |

""" 

def main(argv):
    bcf_in = VariantFile(argv[0])  # auto-detect input format
    print(bcf_in.header, end='')
    for site in bcf_in.fetch():
        keep_site = 0 # default option is to remove SNP
        for sample, rec in site.samples.items():
            if max(rec.get('GP')[1:]) > 0.9:
                keep_site = 1 # do not remove SNP if either het or non-ref homo is greater than .9 for any sample
        if keep_site:
            print(site, end='')
                      
def count_genotypes(fields):
    genotype_counts = [0,0,0,0,0,0]
    for field in fields:
        genotype = field.split(':')[0]
        genotype = genotype.replace("|", "\\")
        if genotype == "0/0":
           genotype_counts[0] += 1
        elif genotype == "0/1" or genotype == "1/0" :
           genotype_counts[1] += 1
        elif genotype == "1/1":
           genotype_counts[2] += 1
        elif genotype == "0/." or genotype == "./0" :
           genotype_counts[3] += 1
        elif genotype == "1/." or genotype == "./1" :
           genotype_counts[4] += 1
        elif genotype == "./.":
           genotype_counts[5] += 1
        else:
            warnings.warn("genotype %s not recognised" % genotype)
    #genotype_counts = [str(i) for i in genotype_counts]
    return genotype_counts          
    

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
           
if __name__ == "__main__":
    global usage
    usage = "python FilterVCF.py VCF_file.vcf.gz | bcftools view -Ob -o VCF_file.filter.vcf.gz (must have corresponding VCF_file.vcf.gz/tbi file)"
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
