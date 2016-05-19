#!/usr/bin/env python

import sys
import subprocess
import fileinput
import warnings
import mysql.connector 
from mysql.connector import Error
from string import maketrans

"""
Parse GTF from Tophat and add to DB
""" 

def main(sample_id):
    try:
       conn = mysql.connector.connect(host='localhost',
                                       database='FetalRNAseq',
                                       #autocommit=True ,
                                       user='root'
                                       )
            
       cursor = conn.cursor()
             
       for line in fileinput.input([]):
           line = line.strip()
           
           parsed = parse_GTF(line.split('\t'))
           keys = parsed.keys()
           if parsed['feature'] == 'transcript':
               try:
                   cursor.execute('INSERT INTO `CufflinksGTF`(frac,full_read_support,cov,conf_lo,FPKM,feature,start,transcript_id,gene_id,conf_hi,seqid,end,score,strand,source,sample_id) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', (
                          parsed['frac'],
                          parsed['full_read_support'],
                          parsed['cov'],
                          parsed['conf_lo'],
                          parsed['FPKM'],
                          parsed['feature'],
                          parsed['start'],
                          parsed['transcript_id'],
                          parsed['gene_id'],
                          parsed['conf_hi'],
                          parsed['seqid'],
                          parsed['end'],
                          parsed['score'],
                          parsed['strand'],
                          parsed['source'],
                          sample_id,))
               except mdb.IntegrityError, e:
                  warnings.warn("%s" % e)
                  pass

           elif parsed['feature'] == 'exon':
               try:
                   cursor.execute('INSERT INTO `CufflinksGTF`(frac,exon_number,cov,conf_lo,FPKM,feature,start,transcript_id,gene_id,conf_hi,seqid,end,score,strand,source,sample_id) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', (
                          parsed['frac'],
                          parsed['exon_number'],
                          parsed['cov'],
                          parsed['conf_lo'],
                          parsed['FPKM'],
                          parsed['feature'],
                          parsed['start'],
                          parsed['transcript_id'],
                          parsed['gene_id'],
                          parsed['conf_hi'],
                          parsed['seqid'],
                          parsed['end'],
                          parsed['score'],
                          parsed['strand'],
                          parsed['source'],
                          sample_id,))
               except mdb.IntegrityError, e:
                  warnings.warn("%s" % e)
                  pass

           else:
               warnings.warn("feature %s not recognised" % parsed['feature'])
           cursor.execute("SELECT * FROM CufflinksGTF")
           results = cursor.fetchall()
    except Error as e:
        print(e)
 
    finally:
        conn.commit()
        cursor.close()
        conn.close()
    
    
    
    
    
def parse_GTF (fields):
  tags = {}
  for attribute in fields[8].split(";")[:-1]:
    attribute = attribute.strip()
    tags[attribute.split(" ")[0]] = " ".join(attribute.split(" ")[1:]).replace('"','')
  try:
    tags['frame'] = int(fields[7])
  except ValueError:
    if fields[7] == '.':
      tags['frame'] = ''
    else:
      sys.exit("frame %s not recognized. Must be 1, 2, 3 or ." % fields[5]) 
  if fields[6] == '-' or fields[6] == 0 or fields[6] == -1:
    tags['strand'] = 0
  elif fields[6] == '+' or fields[6] == 1:
    tags['strand'] = 1
  elif fields[6] == '.':
    tags['strand'] = ''
  else:  
    sys.exit("strand %s not recognized. Must one of +, -, 1, -1 or 0" % fields[6])      
  try:
    tags['score'] = float(fields[5])
  except ValueError:
    if fields[5] == '.':
      tags['score'] = ''
    else:
      sys.exit("score %s not recognized. Must be a number" % fields[5])
  try:
    tags['end'] = int(fields[4])
  except ValueError:
    sys.exit("score %s not recognized. Must be a positive integer" % fields[4])
  try:
    tags['start'] = int(fields[3])
  except ValueError:
    sys.exit("score %s not recognized. Must be a positive integer" % fields[3])
  tags['feature'] = fields[2]
  tags['source'] = fields[1]
  tags['seqid'] = fields[0]
  return tags
    
def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    sample_id = sys.argv[1]
    main(sample_id)
