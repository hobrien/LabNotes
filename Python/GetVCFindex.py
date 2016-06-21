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
             
       for line in fileinput.input():
           (pos, sentrix1, sentrix2) = line.strip().split('_')
           cursor.execute("SELECT BrainBankID FROM PC_analysis WHERE Sentrix_Full = %s", ('_'.join((sentrix1,sentrix2)),))
           rows = cursor.fetchall()
           assert len(rows) == 1
           print '\t'.join((rows[0][0], str(int(pos)+9)))
  
    except Error as e:
        print(e)
 
    finally:
        cursor.close()
        conn.close()
    
    

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)


           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    main(sys.argv[1:])
    
