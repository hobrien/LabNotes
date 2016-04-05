#!/usr/bin/env python

import fileinput

def main():
    for line in fileinput.input():
        result = parse_blast_stats(line)
        if result['evalue'] == 0.0:
            print write_bed12(result)

def write_bed12(result):
    if result['sstart'] > result['send']:
        (result['sstart'], result['send']) = (result['send'], result['sstart'])
    chrom = result['sseqid']
    chromStart = str(result['sstart'] - 1)
    chromEnd = str(result['send'])
    name = '|'.join((chrom, chromStart, chromEnd))
    score = str(result['bitscore'])
    strand = '-'
    if result['strand'] == 1:
        strand = '+'
    thickStart = str(result['sstart'] - 1)
    thickEnd = str(result['send'])
    itemRgb = '0'
    blockCount = '1'
    blockSizes = str(result['send']-result['sstart']+1)
    blockStarts = '0'
    return '\t'.join((chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts))



def parse_blast_stats(row, column_names='qseqid sseqid pident length mismatch gaps start end sstart send evalue bitscore'):
  """This will convert stats from blast hits to the correct numeric format and determine strand
  start and end coordinates are not reversed for negative strand results"""
  row = row.split()
  result = {}
  for field in column_names.split():
    try:
      if field in ('qlen', 'slen', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qframe', 'sstart', 'send', 'sframe', 'sstrand'):
        result[field] = int(row.pop(0))
      elif field in ('pident', 'evalue', 'bitscore'):
        result[field] = float(row.pop(0))
      else:  
        result[field] = row.pop(0)
    except IndexError:
      sys.exit("number of columns does not match specified file format. Please recheck column headers")
  if len(row) > 0:
    sys.exit("number of columns does not match specified file format. Please recheck column headers")
  if result['sstart'] < result['send']:
    result['strand'] = 1
  else:
    result['strand'] = -1
    
  
  return result    
  
if __name__ == "__main__":
   main()
 