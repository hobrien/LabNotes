#!/usr/bin/env python

import sys
import fileinput
import warnings

"""
Parse GTF from Tophat and add to DB
""" 

def main(sample_id):
    print "CREATE TABLE `CufflinksGTF` ("
    print "`id` int(11) unsigned NOT NULL AUTO_INCREMENT,"     
    for line in fileinput.input([]):
         line = line.strip()
           
         parsed = parse_GTF(line.split('\t'))
         keys = parsed.keys()
         for field in keys:
            if field in ('start', 'end', 'strand', 'frame'):
              print "`%s` int(11) DEFAULT NULL," % field
            elif field in ('full_read_support', 'feature', 'transcript_id', 'gene_id', 'seqid', 'source'):
              print "`%s` varchar(255) DEFAULT NULL," % field
            else:
              print "`%s` float DEFAULT NULL," % field
         print "PRIMARY KEY (`id`),"
         for field in ('sample_id', 'transcript_id', 'gene_id', 'seqid', 'start', 'end'):
              print "KEY `%s` (`%s`)," % (field, field)
         print "`sample_id` varchar(255) DEFAULT NULL" 
         print ") ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;"
        # print 'hi ' * (len(keys)+1)
         print "cursor.execute('INSERT INTO `CufflinksGTF`(%s) VALUES(%s)', (" % (','.join(keys + ['sample_id']), '%s,' * len(keys) + '%s')
         for key in keys:
             print "                          parsed['%s']," % key
         print "                          sample_id,)"    
         break
         
def parse_GTF (fields):
  tags = {}
  for attribute in fields[8].split(";")[:-1]:
    attribute = attribute.strip()
    tags[attribute.split(" ")[0]] = " ".join(attribute.split(" ")[1:]).replace('"','')
  try:
    tags['frame'] = int(fields[7])
  except ValueError:
    if fields[7] == '.':
      tags['frame'] = 'NULL'
    else:
      sys.exit("frame %s not recognized. Must be 1, 2, 3 or ." % fields[5]) 
  if fields[6] == '-' or fields[6] == 0 or fields[6] == -1:
    tags['strand'] = 0
  elif fields[6] == '+' or fields[6] == 1:
    tags['strand'] = 1
  elif fields[6] == '.':
    tags['strand'] = 'NULL'
  else:  
    sys.exit("strand %s not recognized. Must one of +, -, 1, -1 or 0" % fields[6])      
  try:
    tags['score'] = float(fields[5])
  except ValueError:
    if fields[5] == '.':
      tags['score'] = 'NULL'
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
