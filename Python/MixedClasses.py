#!/usr/bin/env python

import sys
import subprocess
import fileinput
import warnings
import mysql.connector 
from mysql.connector import Error


""" for all cuffcompare results with mixed (.) classes, make count of each class in 
corresponding tmap files
this will involve using the following SQL command to fetch the class for each column
after the first 4:

SELECT Cufflinks.class FROM Cufflinks, QueryNumber WHERE Cufflinks.sample_id = 
   QueryNumber.sample_id AND QueryNumber.id = %s AND Cufflinks.transcript_id = %s

entries look like "q24:CUFF.20|CUFF.20.1|100|0.226699|0.131164|0.322235|3.561109|1562"
so id = entry.split(':')[0][1:] and transcript_id = entry.split('|')[1]
""" 

def main(args):
    try:
       conn = mysql.connector.connect(host='localhost',
                                       database='FetalRNAseq',
                                       user='root')
            
       cursor = conn.cursor()
             
       matchingSNPs = 0
       for line in fileinput.input():
           row = line.strip().split('\t')
           counts = {'=':0, 'i':0, 'r':0, 'u':0, 'j':0, 'x':0, 'e':0, 'p':0, 'o':0, 'c':0, 's':0}
           if row[3] != '.':
               continue
           for entry in row[4:]:
               if entry == '-':
                   continue
               id = entry.split(':')[0][1:]
               transcript_id = entry.split('|')[1]
               cursor.execute("SELECT Cufflinks.class FROM Cufflinks, QueryNumber WHERE Cufflinks.sample_id = QueryNumber.sample_id AND QueryNumber.id = %s AND Cufflinks.transcript_id = %s", (id, transcript_id,))
               rows = cursor.fetchall()
               if len(rows) == 0:
                   warnings.warn("no entry for q%s|%s" % (id, transcript_id))
                   continue
               elif len(rows) > 1:
                   warnings.warn("Multiple entires for q%s|%s" % (id, transcript_id))
               counts[rows[0][0]] += 1
           print '\t'.join(str(i) for i in (row[0], counts['='], counts['i'], counts['r'], counts['u'], counts['j'], counts['x'], counts['e'], counts['p'], counts['o'], counts['c'], counts['s']))
  
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
    
