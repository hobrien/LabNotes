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
               coords = GetCoords(cursor, id)
               print "%s:%i-%i" % (coords[0], coords[1] + 1, coords[2])
  
    except Error as e:
        print(e)
 
    finally:
        cursor.close()
        conn.close()
    
    
def GetCoords(cursor, id):
    cursor.execute("SELECT chrom, chromStart, chromEnd FROM snp144 WHERE name = %s", (id,))
    rows = cursor.fetchall()
    ref_row = 0
    try:
        assert len(rows) == 1
    except AssertionError:
        if len(rows) == 0:
            warnings.warn("No SNP named %s in DB" % (id))
        else:
            ref_hap = []  # DB in includes positions of SNPs on alternate haplotypes (eg chr6_apd_hap1). Make a list of rows that have chrom values like 'chr6'
            for index in range(len(rows)):
                if len(rows[index][0].split('_')) == 1:
                    ref_hap.append(index)
            try:
                assert len(ref_hap) == 1
            except AssertionError:
                if len(ref_hap) == 0:
                    warnings.warn("No SNP on reference haplotype named %s in DB" % (id))
                else:            
                    warnings.warn("%i rows on reference haplotype matching name %s" % (len(rows), id))
                sys.exit()    
            ref_row = ref_hap[0]        
    return  rows[ref_row]

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)


           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
    
