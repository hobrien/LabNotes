#!/usr/bin/env python

import sys
import subprocess
import fileinput
import warnings

#Read count in column 4 includes skips caused by spliced alignments. I need to screen these out and report only ref matches and mismatches
def main(args):
       for line in fileinput.input():
           fields = line.split()
           try:
               bases = fields[4]
               ref = bases.count(',')+bases.count('.')
               alt = int(fields[3]) - ref - bases.count('<') - bases.count('>')
               print '\t'.join(fields[:3]+ [str(ref), str(alt)])    
           except IndexError: # no pileup data
               print '\t'.join(fields[:3]+ ['0', '0'])    
  
def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)


           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
    
