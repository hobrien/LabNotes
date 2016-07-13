#!/usr/bin/env python

import sys
import fileinput
import warnings
from string import maketrans

"""

I need to read a line from a vcf and output the genotypes with the BrainBank IDs
should look like this:
chr2	184936178	A|A	15240
chr2	184936178	T|A	15533
...

column 1 and 2 are columns 1 and 2 of the vcf
The next column can be derived by substituting column 4 for 0 and column 5 for 1 in field.split(':')[0]

intab='01'
outtab=fields[3]+fields[4]
trantab = maketrans(intab, outtab)
genotype = field.split(':')[0].translate(trantab)

This is run on fields 10 + 
The IDs are culled from ~/LabNotes/VCFindex.txt (need to subtract 1 for zero-based)

        
""" 

def main(argv):
    with open(argv[0], 'r') as index_file:
        index_file = open(argv[0], "rw+")
        VCF_index = index_file.readlines()
        for line in fileinput.input([]):
           line = line.strip()
           get_genotypes(VCF_index, line)

def get_genotypes(VCF_index, VCF_line):
    fields = VCF_line.split()
    intab='01'
    try:
        outtab=fields[3]+fields[4]
    except IndexError:
        warnings.warn("Not enough columns in VCF. This is usually because the headeris included. Be sure to use -H option for bcftools")
        sys.exit("Usage: %s" % usage)
    trantab = maketrans(intab, outtab)
    for sample in VCF_index:
        (sampleID, index) = sample.split('\t')
        index = int(index) - 1 # one-based to zero based
        try:
            genotype = fields[index]
        except IndexError:
            sys.exit("the number of samples in the index file does not match the number of fields in the VCF")
        genotype = genotype.split(':')[0].translate(trantab)
        print '\t'.join((fields[0], fields[1], genotype, sampleID))
    
def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
           
if __name__ == "__main__":
    global usage = "bcftools view -H -r chrX:XXX XXX.vcf | python GetGenotypes.py VCF_index.txt"
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
