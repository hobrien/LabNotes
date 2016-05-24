#!/usr/bin/env python

import sys
import fileinput
import warnings

"""
Parse GTF from Tophat and add to DB
""" 

def main(argv):
       id = 0
       fixed_fields = ('seqid', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame')
       #I'm going to write the features to aa file with an arbitrary name and the other fields to STDOUT
       with open("features.csv", 'w') as feature_fh:
         for line in fileinput.input():
           line = line.strip()
           parsed = parse_GTF(line.split('\t'))
           if not parsed:
               continue
           keys = parsed.keys()
           id += 1
           print ','.join([str(id)] + [str(parsed[x]) for x in fixed_fields])
           for key in keys:
               if key not in fixed_fields:
                   feature_fh.write(','.join([str(x) for x in (id, key, parsed[key])])+'\n')

           
    
    
    
    
    
def parse_GTF (fields):
  tags = {}
  if fields[0][:2] == '##':
      return ''
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
    main(sys.argv[1:])
