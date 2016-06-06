#!/usr/bin/env python

import sys
import subprocess
import fileinput
import warnings
import mysql.connector 
from mysql.connector import Error


""" for each Gene symbol in input, get homolog ID from DB, check that there are a total of
2 entries for that homologID, then replace thr Gene symbol with the symbol entry from the 
DB for the human homolog

These are the SQL commands needed (could probably be combined, but I don't know how):
SELECT MouseHumanHomo.`HomoloGene ID` FROM MouseHumanHomo, MGI WHERE MouseHumanHomo.`Mouse MGI ID` = MGI.`MGI Gene ID` AND MGI.Input = %s;

SELECT COUNT(*) FROM MouseHumanHomo WHERE `HomoloGene ID` = %s;

SELECT Symbol FROM MouseHumanHomo WHERE `HomoloGene ID` = 9170 AND `NCBI Taxon ID` = %s;

""" 

def main():
    try:
       conn = mysql.connector.connect(host='localhost',
                                       database='FetalRNAseq',
                                       user='root')
            
       cursor = conn.cursor()
             
       for line in fileinput.input([]):
           line = line.strip()
           gene_id = line.split('\t')[0]
           if gene_id == 'Gene symbol':
               print '\t'.join([gene_id] + line.split('\t')[2:])
               continue
           cursor.execute("SELECT `MGI Gene ID` FROM MGI WHERE MGI.Input = %s AND (MGI.`Input type` = 'current symbol' OR MGI.`Input type` = 'old symbol' OR MGI.`Input type` = 'related synonym' OR MGI.`Input type` = 'GenBank' OR MGI.`Input type` = 'synonym')", (gene_id,))
           #print "SELECT Cufflinks.class FROM Cufflinks WHERE Cufflinks.sample_id = %s AND Cufflinks.transcript_id = %s" % (sample_id, transcript_id,)
           rows = cursor.fetchall()
           if len(rows) == 0:
                   warnings.warn("no MGI entry for %s" % (gene_id, ))
                   continue
           elif len(rows) > 1:
               for row in rows[1:]:
                   if not row[0] == rows[0][0]: 
                       warnings.warn("Multiple conflicting MGI entires for %s" % (gene_id, ))
                       continue
           MGI_id = rows[0][0]
 
           cursor.execute("SELECT `HomoloGene ID` FROM MouseHumanHomo WHERE `Mouse MGI ID` = %s", (MGI_id,))
           #print "SELECT Cufflinks.class FROM Cufflinks WHERE Cufflinks.sample_id = %s AND Cufflinks.transcript_id = %s" % (sample_id, transcript_id,)
           rows = cursor.fetchall()
           if len(rows) == 0:
                   #warnings.warn("no Homolog entry for %s" % (gene_id, ))
                   continue
           elif len(rows) > 1:
                   #warnings.warn("Multiple Homolog entires for %s" % (gene_id, ))
                   continue           
           homo_id = rows[0][0]

           cursor.execute("SELECT COUNT(*) FROM MouseHumanHomo WHERE `HomoloGene ID` = %s", (homo_id,))
           if cursor.fetchone()[0] != 2:
               continue

           cursor.execute("SELECT Symbol FROM MouseHumanHomo WHERE `HomoloGene ID` = %s AND `NCBI Taxon ID` = 9606", (homo_id,))
           rows = cursor.fetchall()
           if len(rows) != 1:
               continue
           print '\t'.join([rows[0][0]] + line.split('\t')[2:])
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
    
