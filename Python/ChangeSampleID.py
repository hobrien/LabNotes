#!/usr/bin/env python

import sys
import subprocess
import mysql.connector 
from mysql.connector import Error

# read header using bcftools, lookup sample IDs in DB, then write new header with BrainIDs
 

def main(args):
    vcf = args[0]
    with open('temp.head', 'w') as temp_header:
        proc = subprocess.Popen(['bcftools', 'view', '-h',  vcf],stdout=subprocess.PIPE)
        for line in iter(proc.stdout.readline,''):
            line = line.strip()
            if line[:6] == '#CHROM':
               line = ConvertIDs(line)
            temp_header.write(line + '\n')


def ConvertIDs(line):
    line = line.strip()
    headers = line.split('\t')
    for sampleID in range(9, len(headers)):
        headers[sampleID] = BrainID(headers[sampleID][3:]) 
    line ='\t'.join(headers)
    return line
    
def BrainID(Sentrix):
    try:
       conn = mysql.connector.connect(host='localhost',
                                       database='FetalRNAseq',
                                       user='root')
            
       cursor = conn.cursor()
       cursor.execute("SELECT BrainBankID FROM PC_analysis WHERE Sentrix_Full = %s", (Sentrix,))
       rows = cursor.fetchall()
       try:
           assert cursor.rowcount == 1
       except AssertionError:
           if cursor.rowcount == 0:
               sys.exit("%s not recognised" % Sentrix)
           else:
               sys.exit("%i rows matching %s" % (cursor.rowcount, Sentrix))
       return  rows[0][0]
             
    except Error as e:
        print(e)
 
    finally:
        cursor.close()
        conn.close()
        

           
if __name__ == "__main__":
    main(sys.argv[1:])
    
