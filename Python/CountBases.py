#!/usr/bin/env python

import sys
import subprocess
import fileinput
import warnings
import mysql.connector 
from mysql.connector import Error

# read header using bcftools, lookup sample IDs in DB, then write new header with BrainIDs
 

def main(args):
       for line in fileinput.input():
           fields = line.split()
           bases = fields[4]
           counts = bases.count(',')+bases.count('.')
           print '\t'.join(fields[:4]+ [str(counts)])    
               
  
def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)


           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
    
