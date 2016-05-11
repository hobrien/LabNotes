#!/usr/bin/env python

import sys
import subprocess
import fileinput
import warnings
import mysql.connector 
from mysql.connector import Error


""" for each entry in GTF files, find corresponding class code from DB of tracking file
and add it as an atribute at the end of the line

This is the SQL command needed:
SELECT Cufflinks.class FROM Cufflinks WHERE Cufflinks.id = %s AND Cufflinks.transcript_id = %s

transcript_id looks to be the 12th column so this should work (also need to remove quotes and semicolon)
transcript_id = line.split()[11][1:-2]

It is probably going to be easiest to just feed the sample_id from the command line

This should work for the output:

print '%s; class "%s"' % (line, class)
""" 

def main(sample_id):
    try:
       conn = mysql.connector.connect(host='localhost',
                                       database='FetalRNAseq',
                                       user='root')
            
       cursor = conn.cursor()
             
       for line in fileinput.input([]):
           line = line.strip()
           transcript_id = line.split()[11][1:-2]
           cursor.execute("SELECT Cufflinks.class FROM Cufflinks WHERE Cufflinks.sample_id = %s AND Cufflinks.transcript_id = %s", (sample_id, transcript_id,))
           #print "SELECT Cufflinks.class FROM Cufflinks WHERE Cufflinks.sample_id = %s AND Cufflinks.transcript_id = %s" % (sample_id, transcript_id,)
           rows = cursor.fetchall()
           if len(rows) == 0:
                   #warnings.warn("no entry for q%s|%s" % (sample_id, transcript_id))
                   print line
                   continue
           elif len(rows) > 1:
                   warnings.warn("Multiple entires for q%s|%s" % (sample_id, transcript_id))
                   print line
                   continue
           transfrag_class = rows[0][0]
           print '%s class "%s";' % (line, transfrag_class)
    except Error as e:
        print(e)
 
    finally:
        cursor.close()
        conn.close()
    
    
def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    sample_id = sys.argv[1]
    main(sample_id)
    
