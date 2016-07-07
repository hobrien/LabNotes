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
       index = 0      
       for line in fileinput.input():
           (pos, sentrix1, sentrix2) = line.strip().split('_')
           index += 1
           sid = '_'.join((sentrix1,sentrix2))
           cursor.execute("SELECT BrainBankID FROM PC_analysis WHERE Sentrix_Full = %s", (sid,))
           rows = cursor.fetchall()
           if len(rows) == 0:
               warnings.warn("No DB entry for %s" % sid)
               sys.exit()
           elif len(rows) > 1:    
               warnings.warn("Multiple DB entry for %s" % sid)
               sys.exit()
           print '\t'.join((rows[0][0], str(int(index)+9)))
  
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
    
