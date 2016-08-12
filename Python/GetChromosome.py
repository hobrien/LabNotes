#!/usr/bin/env python

import sys
import subprocess
import fileinput
import warnings
import mysql.connector 
from mysql.connector import Error

# Read a list of geneIDs and add chromosome
 

def main():
    try:
       conn = mysql.connector.connect(host='localhost',
                                       database='FetalRNAseq',
                                       user='root')
            
       cursor = conn.cursor()
             
       header = 0
       for line in fileinput.input():
           fields = line.strip().split('\t')
           if header == 0:
               print '\t'.join(['Chromosome', 'Gene'] + fields)
               header = 1
           else:    
               (chrom, gene) = GetChr(cursor, fields[0])
               print '\t'.join([chrom, gene] + fields)
  
    except Error as e:
        print(e)
 
    finally:
        cursor.close()
        conn.close()
    
    
def GetChr(cursor, gene_id):
    cursor.execute("SELECT DISTINCT GencodeGTF.seqid, GencodeFeatures.value FROM GencodeFeatures, GencodeGTF WHERE GencodeGTF.id = GencodeFeatures.id AND GencodeFeatures.feature = 'gene_name' AND GencodeFeatures.id IN (SELECT GencodeGTF.id FROM GencodeGTF, GencodeFeatures WHERE GencodeGTF.id = GencodeFeatures.id AND GencodeFeatures.Value = %s)", (gene_id,))
    rows = cursor.fetchall()
    ref_row = 0
    try:
        assert len(rows) == 1
    except AssertionError:
        if len(rows) == 0:
            warnings.warn("No feature with geneID %s in DB" % (gene_id))
            rows = [['NA']]
        else:
            warnings.warn("Feature with geneID %s on multiple chromosomes" % (gene_id))           
    return  rows[0]

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)


           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    main()
    
