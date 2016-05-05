#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
#path=${1%/*}
#sampleID=${path##*/}

echo "Starting Cuffcompare"
cp /c8000xd3/rnaseq-heath/Cufflinks/15468/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/15468.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/15533/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/15533.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/15768/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/15768.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/16286/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16286.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/16438/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16438.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/16488/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16488.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/16840/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16840.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/16929/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16929.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/16972/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16972.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/17049/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17049.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/17054/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17054.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/17068/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17068.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/17081/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17081.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/17087/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17087.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/17109/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17109.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/17115/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17115.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/17130/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17130.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/17701/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17701.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/17812/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17812.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/18294/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18294.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/18349/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18349.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/18596/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18596.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/18655/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18655.gtf
cp /c8000xd3/rnaseq-heath/Cufflinks/18687/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18687.gtf
cd /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/
cuffcompare -V -o Combined \
  -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
  -s /c8000xd3/rnaseq-heath/Ref/chroms/ \
  15468.gtf 15533.gtf 15768.gtf 16286.gtf 16438.gtf 16488.gtf 16840.gtf 16929.gtf 16972.gtf \
  17049.gtf 17054.gtf 17068.gtf 17081.gtf 17087.gtf 17109.gtf 17115.gtf 17130.gtf 17701.gtf \
  17812.gtf 18294.gtf 18349.gtf 18596.gtf 18655.gtf 18687.gtf
echo "Finished Cuffcompare"
