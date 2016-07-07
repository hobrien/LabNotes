#!/usr/bin/env python

import sys
import subprocess
import fileinput
import warnings
import mysql.connector 
from mysql.connector import Error


"""The file that Nick sent me with the Sentrix IDs only included one batch of 96
I've found another file with all of the data, so I need to update what is in the DB
with the missing info (only the sentrix IDs and BrainBank IDs). I will also do a check to
make sure that BrainBank IDs that are already in the DB have matching Sentrix numbers
""" 

def main():
    try:
       conn = mysql.connector.connect(host='localhost',
                                      database='FetalRNAseq',
                                      autocommit=True ,
                                      user='root')
            
       cursor = conn.cursor()
             
       for line in fileinput.input([]):
           line = line.strip()
           (bbid, sid) = line.split('\t')
           cursor.execute("SELECT Sentrix_Full FROM PC_analysis WHERE BrainBankID = %s", (bbid,))
           rows = cursor.fetchall()
           if len(rows) == 0:
               print "INSERT INTO PC_analysis (BrainBankID, Sentrix_Full) VALUES (%s, %s)" % (bbid,sid)
               cursor.execute("INSERT INTO PC_analysis (BrainBankID, Sentrix_Full) VALUES (%s, %s)", (bbid,sid,))
           elif len(rows) == 1:
               if rows[0][0] != sid:
                   warnings.warn("DB entry for %s (%s) does not match input (%s)" % (bbid, rows[0][0], sid))
           else:
                   warnings.warn("Multiple entires for q%s|%s" % (sample_id, transcript_id))
    except Error as e:
        print(e)
 
    finally:
        cursor.close()
        conn.close()
    
    
def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    main()
    
