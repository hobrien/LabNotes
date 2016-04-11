#!/usr/bin/env python

import sys
import subprocess
import fileinput
import warnings
import mysql.connector 
from mysql.connector import Error

# read header using bcftools, lookup sample IDs in DB, then write new header with BrainIDs
 

def main(args):
    try:
       conn = mysql.connector.connect(host='localhost',
                                       database='FetalRNAseq',
                                       user='root')
            
       cursor = conn.cursor()
             
       matchingSNPs = 0
       for line in fileinput.input():
           for id in line.split():
               print '\t'.join(str(x) for x in GetCoords(cursor, id))
  
    except Error as e:
        print(e)
 
    finally:
        cursor.close()
        conn.close()
    
    
def GetCoords(cursor, id):
    cursor.execute("SELECT chrom, chromStart, chromEnd FROM snp144 WHERE name = %s", (id,))
    rows = cursor.fetchall()
    try:
        assert len(rows) == 1
    except AssertionError:
        if len(rows) == 0:
            warnings.warn("No SNP at %s position %s" % (chromosome, position))
        else:
            warnings.warn("%i rows matching %s position %s" % (len(rsID), chromosome, position))
        sys.exit()    
    return  rows[0]

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)


           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
    
