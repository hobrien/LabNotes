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
           (chromosome, position, id) = line.strip().split('\t')
           if id != GetID(cursor, chromosome, position):
               warnings.warn("id %s at %s position %s does not match DB id (%s)" % (id, chromosome, position, GetID(cursor, chromosome, position)))
           matchingSNPs += 1
           if matchingSNPs % 1000 == 0:
               warnings.warn("matched %i SNPs" % matchingSNPs)
  
    except Error as e:
        print(e)
 
    finally:
        cursor.close()
        conn.close()
    

def ConvertIDs(line):
    line = line.strip()
    headers = line.split('\t')
    for sampleID in range(9, len(headers)):
        headers[sampleID] = BrainID(headers[sampleID][3:]) 
    line ='\t'.join(headers)
    return line
    
def GetID(cursor, chromosome, position):
    if chromosome[:3] != 'Chr':
        chromosome = "Chr%s" % chromosome
    cursor.execute("SELECT ID FROM dbSNP146 WHERE Chr = %s AND end = %s", (chromosome, position,))
    rows = cursor.fetchall()
    try:
        assert cursor.rowcount == 1
    except AssertionError:
        if cursor.rowcount == 0:
            warnings.warn("No SNP at %s position %s" % (chromosome, position))
            rows = [['NA']]
        else:
            warnings.warn("%i rows matching %s position %s" % (cursor.rowcount, chromosome, position))
    return  rows[0][0]

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)


           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
    
