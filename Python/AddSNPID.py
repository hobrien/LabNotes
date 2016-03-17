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
           (chromosome, position, id, ref, alt) = line.strip().split('\t')
           if id != GetID(cursor, chromosome, position, ref, alt):
               warnings.warn("id %s at %s position %s does not match DB id (%s)" % (id, chromosome, position, GetID(cursor, chromosome, position, ref, alt)))
           matchingSNPs += 1
           if matchingSNPs % 1000 == 0:
               warnings.warn("matched %i SNPs" % matchingSNPs)
  
    except Error as e:
        print(e)
 
    finally:
        cursor.close()
        conn.close()
    
    
def GetID(cursor, chromosome, position, ref, alt):
    if chromosome[:3] != 'Chr':
        chromosome = "Chr%s" % chromosome
    cursor.execute("SELECT name, observed FROM snp144 WHERE chrom = %s AND chromStart = %s", (chromosome, position - 1,))
    rows = cursor.fetchall()
    rsID = []
    for row in rows:
        db = row[1].encode('ascii').replace('-', '.').split('/')
        db.sort()
        vcf = [ref, alt]
        vcf.sort()
        if db == vcf:
            rsID.append(row[0])
    try:
        assert len(rsID) == 1
    except AssertionError:
        if len(rsID) == 0:
            warnings.warn("No SNP at %s position %s" % (chromosome, position))
            rsID = ['NA']
        else:
            warnings.warn("%i rows matching %s position %s" % (len(rsID), chromosome, position))
    return  rsID[0]

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)


           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
    
