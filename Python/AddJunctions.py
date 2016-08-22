#!/usr/bin/env python

import sys
import subprocess
import fileinput
import warnings
import mysql.connector 
from mysql.connector import Error
from string import maketrans

"""
BED columns are as follows:
   chrom        The name of the chromosome.
   chromStart   The starting position of the junction-anchor. This includes the maximum overhang for the junction on the left. For the exact junction start add blockSizes[0].
   chromEnd     The ending position of the junction-anchor. This includes the maximum overhang for the juncion on the left. For the exact junction end subtract blockSizes[1].
   name         The name of the junctions, the junctions are labeled 'partial_novel' or 'complete_novel'
   score        The number of reads supporting the junction.
   strand       Defines the strand - either '+' or '-'. This is set to '.' for junction_annotation.py output
   thickStart   Same as chromStart.
   thickEnd     Same as chromEnd.
   itemRgb      RGB value 
   blockCount   The number of blocks, 2 by default.
   blockSizes   A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
   blockStarts  A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
    
DB fields
    id         auto-incremented novel intiger
    junc_id    novel id for each junction position
    sample_id  id of sample (junct_id + sample_id should be unique)
    start      chromStart
    end        chromEnd
    strand     strand
    chr        chrom
    status     name
    score      score
    int_start  intron start (chromStart + blockSizes[0] + 1
    int_end    intron_end (chromEnd - blockSizes[1] = chromStart + blockStarts[1]
    position   key on chr, int_Start and int_end
    junction   unique key on junc_id and sample_id (will raise error if multiple junctions with the same id from the same sample
    

- The output from annotate_junctions.py does not give the correct block sizes (10 bp instead of the maximum length of the read)
- It also does not give strand info
- I thought about using the junctions.bed file from TopHat (which does have all the correct info) and identifying novel junctions by searching the
    positions in the reference GTF, but this would miss exon skipping events
- I also thought about using the output form Cufflinks with the 'j' class, but I think this would only catch novel junctions in known genes
- I think cufflinks might also classify junction read-though events as splice variants ('j')
- I only want variants supported by spliced reads
- It should be easy enough to match up the junctions from annotate_junctions.py and update the block size (and chromStart/chromEnd + blockStarts) and strand fields if I need to
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


def parse_BED (line):
  fields = line.split('\t')
  if len(fields) != 12:
      sys.exit("this is not a BED12 file. There are %s columns" % len(fields))   
  if fields[9] != 2:
      sys.exit("this is not a junctions file. There should be 2 blocks, but there are %s" % fields[9])          
  blockSizes = fields[10].split(',')
  if len(blockSizes) != 2:
      sys.exit("this is not a properly formatted BED file. There are 2 blocks but %s blockSizes" % len(blockSizes))
  blockStarts = fields[11].split(',')
  if len(blockStarts) != 2:
      sys.exit("this is not a properly formatted BED file. There are 2 blocks but %s blockStarts" % len(blockSizes))
  if blockStarts[0] != 0:
      sys.exit("this is not a properly formatted BED file. The first block starts %s bases from the chromStart" % blockStarts[0])
  if fields[2] != int(fields[1]) + int(blockStarts[1]) + int(blockSizes[1])
      sys.exit("this is not a properly formatted BED file. The second block ends at position %i but chromEnd is %i" % (int(fields[1]) + int(blockStarts[1]) + int(blockSizes[1]), int(fields[2]))
  
  tags = {}
  tags['chr'] = fields[0]
  tags['start'] = fields[1]
  tags['end'] = fields[2]
  tags['status'] = fields[3]
  tags['score'] = fields[4]
  tags['strand'] = fields[5]
  tags['int_start'] = int(chromStart) + int(blockSizes[0] + 1
  tags['int_end'] = int(chromEnd) - int(blockSizes[1])
  return tags


def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
 
           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    sample_id = sys.argv[1]
    main(sample_id)
