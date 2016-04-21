#!/usr/bin/env python

import sys
import subprocess
import fileinput
import warnings
import mysql.connector 
from mysql.connector import Error
from string import maketrans

# read header using bcftools, lookup sample IDs in DB, then write new header with BrainIDs
 

def main(args):
    try:
       conn = mysql.connector.connect(host='localhost',
                                       database='FetalRNAseq',
                                       user='root')
            
       cursor = conn.cursor()
             
       matchingSNPs = 0
       for line in fileinput.input():
           if line[0] == '#':
               print line.strip()
               continue
           row = line.strip().split('\t')
           (chromosome, position, vcf_id, ref, alt) =  row[:5]
           if alt == '.':
               continue
           db_id = GetID(cursor, chromosome, position, ref, alt, id)
           if db_id != 'NA':
               row[2] = db_id
           print '\t'.join(row)
           matchingSNPs += 1
           #if matchingSNPs % 1000 == 0:
           #    warnings.warn("matched %i SNPs" % matchingSNPs)
  
    except Error as e:
        print(e)
 
    finally:
        cursor.close()
        conn.close()
    
    
def GetID(cursor, chromosome, position, ref, alt, id):
    if chromosome[:3] != 'Chr':
        chromosome = "Chr%s" % chromosome
    cursor.execute("SELECT name, observed, exceptions FROM snp144 WHERE chrom = %s AND chromStart = %s AND chromEnd = %s", (chromosome, int(position)-1, int(position)-1+len(ref),))
    rows = cursor.fetchall()
    
    rsID = []
    if len(rows) == 0:
        return 'NA'
    elif len(rows) == 1:
        rsID = [rows[0][0]]
    else:    
        for row in rows:
            db = row[1].encode('ascii').replace('-', '.').split('/')
            if ref in db:
                if alt in db:
                    rsID.append(row[0])
            elif RevCom(ref) in db and RevCom(alt) in db:
                rsID.append(row[0])
    try:
        assert len(rsID) == 1
    except AssertionError:
        if len(rsID) == 0:
            warnings.warn("None of the %i SNPs at %s position %s match %s/%s" % (len(rows), chromosome, position, ref, alt))
            rsID = ['NA']
        elif 'DuplicateObserved' not in rows[0][2]:
            warnings.warn("%i rows matching %s position %s" % (len(rsID), chromosome, position))
    return  rsID[0]

def RevCom(vcf):
    intab = "ACGT"
    outtab = "TGCA"
    trantab = maketrans(intab, outtab)
    vcf = vcf.translate(trantab)[::-1]
    return vcf
    
def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
    
