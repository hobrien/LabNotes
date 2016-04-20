import sys, warnings
#all parts of the names should match up except for the 'R?[12](_001)?' portion before '.fastq.gz'

seq1 = sys.argv[1]
seq2 = sys.argv[2]

def TruncName(seq):
    if '_R1_001.fastq.gz' in seq:
        seq = seq.replace('_R1_001.fastq.gz' , '')
    elif '_R2_001.fastq.gz' in seq:
        seq = seq.replace('_R2_001.fastq.gz' , '')  
    elif '_1.sanfastq.gz' in seq:
        seq = seq.replace('_1.sanfastq.gz' , '')  
    elif '_2.sanfastq.gz' in seq:
        seq = seq.replace('_2.sanfastq.gz' , '')  
    else:
        sys.exit('name %s not recognised' % seq)
    return seq

if TruncName(seq1) != TruncName(seq2):
   warnings.warn("%s and %s do not match!" % (seq1, seq2))
