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
        VCFindex = index_file.readlines()
        for line in fileinput.input():
           line = line.strip()
           print get_genotypes(VCF_index, line)

def get_genotype(VCF_index, VCF_line):
    fields = VCF_line.split('\t')
    for line in VCFindex:
        (sampleID, index) = line.split('\t')
        index -= 1 # one-based to zero based
        genotype = fields[index].split(':')[0].translate(trantab)
        try:
            return '\t'.join((fields[0], fields[1], genotype, sampleID))
        except IndexError:
            sys.exit("the number of samples in the index file does not match the number of fields in the VCF")
    
def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
