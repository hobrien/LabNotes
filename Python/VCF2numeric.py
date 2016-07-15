#!/usr/bin/env python

import sys
import fileinput
import warnings
from string import maketrans

"""
"For each gene-SNP pair, with the SNP encoded by 0,1 and 2 according to the frequency of the minor allele, the association between gene expression g and genotype s is assumed to be linear"

I interpret this to mean that the het is encoded 1 and the minor allele homozygote is encoded 2 (ie; the number of copies of the minor allele)
I assume that I should use the minor allele of the sample, not the global minor allele (this certainly makes it easier)
The VCF is coded according to the reference allele, so I am going to have to count the number of 0s and 1s to figure out which is the minor allele and recode accordingly.


""" 

def main(argv):
    for line in fileinput.input([]):
       line = line.strip()
       try:
           if line[:2] == '##': # skip header lines
               continue
           elif line[0] == '#': # column names
               fields = line.split('\t')
               print '\t'.join((fields[2]) + fields[9:])
        except IndexError:
            continue
        fields = line.split('t')
        
def count_alleles(line):
    

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
           
if __name__ == "__main__":
    global usage
    usage = "bcftools view -H -r chrX:XXX XXX.vcf | python GetGenotypes.py VCF_index.txt"
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
