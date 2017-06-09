#!/usr/bin/env python

import sys
import fileinput
import warnings
from string import maketrans
import re

"""
"For each gene-SNP pair, with the SNP encoded by 0,1 and 2 according to the frequency of the minor allele, the association between gene expression g and genotype s is assumed to be linear"

I interpret this to mean that the het is encoded 1 and the minor allele homozygote is encoded 2 (ie; the number of copies of the minor allele)
I assume that I should use the minor allele of the sample, not the global minor allele (this certainly makes it easier)
The VCF is coded according to the reference allele, so I am going to have to count the number of 0s and 1s to figure out which is the minor allele and recode accordingly.

I'm not sure what to do about missing data. Currently I'm encoding any SNP with a missing allele as homozygous for the major allele.
""" 

def main():
    for line in fileinput.input([]):
       line = line.strip()
       try:
           if line[:2] == '##': # skip header lines
               continue
           elif line[0] == '#': # column names
               fields = line.split('\t')
               print '\t'.join(['id'] + fields[9:])
               continue
       except IndexError:
           continue
       line=line.replace('/', '|')  
       fields = line.split('\t')
       output = [fields[2]]
       if minor_allele(line) == 0:
           minor_homo = '0|0'
           major_homo = '1|1'
       else:
           minor_homo = '1|1'
           major_homo = '0|0'
           
       for field in fields[9:]:
           if field[:3] == major_homo or field[:3].find('.') >= 0:
               output += ['0']
           elif field[:3] == '0|1' or field[:3] == '1|0':
               output += ['1']
           elif field[:3] == minor_homo:
               output += ['2']
           else:
               print field
               sys.exit(field)
       print '\t'.join(output)                  
               
        
def minor_allele(line):
    if len(re.findall(r'(?=(0\|)|(\|0))', line)) <= len(re.findall(r'(?=(1\|)|(\|1))', line)):
        return 0
    else:
        return 1
    

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
           
if __name__ == "__main__":
    global usage
    usage = "bcftools view -H -r chrX:XXX XXX.vcf | python GetGenotypes.py VCF_index.txt"
    warnings.formatwarning = warning_on_one_line
    main()
