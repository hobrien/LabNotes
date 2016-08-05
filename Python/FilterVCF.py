#!/usr/bin/env python

import sys
import fileinput
import warnings
from string import maketrans

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
    print '\t'.join(('Position', '0/0', '0/1', '1/1', '0/.', '1/.', './.'))
    for line in fileinput.input([]):
       line = line.strip()
       try:
           if line[0] == '#': # skip header lines and column names
               continue
       except IndexError:
            continue
       fields = line.split('\t')
       genotype_counts = count_genotypes(fields[9:])
       if sum(genotype_counts[3:]) == 0:
           print line
        
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
    usage = "bcftools view -H -r chrX:XXX XXX.vcf | python GetGenotypes.py VCF_index.txt"
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
