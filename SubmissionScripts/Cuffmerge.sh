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

echo "Starting Cuffmerge"
#cp /c8000xd3/rnaseq-heath/Cufflinks/15468/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/15468.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/15533/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/15533.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/15768/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/15768.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16286/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16286.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16438/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16438.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16488/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16488.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16840/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16840.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16929/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16929.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16972/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16972.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17049/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17049.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17054/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17054.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17068/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17068.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17081/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17081.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17087/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17087.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17109/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17109.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17115/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17115.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17130/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17130.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17701/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17701.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17812/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17812.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18294/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18294.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18349/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18349.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18596/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18596.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18655/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18655.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18687/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18687.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/15641/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/15641.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16024/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16024.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16115/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16115.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16385/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16385.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16428/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16428.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16491/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16491.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16548/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16548.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16810/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16810.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16826/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16826.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17048/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17048.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17053/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17053.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17071/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17071.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17921-l1/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17921-l1.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/15655/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/15655.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16483/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16483.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16640/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16640.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17013/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17013.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17198/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17198.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17229/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17229.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17333/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17333.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17475/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17475.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17543/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17543.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17629/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17629.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17835/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17835.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18249/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18249.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18372/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18372.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18666/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18666.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18983/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18983.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/19043/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/19043.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/19052/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/19052.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/A138/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/A138.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/A226/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/A226.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/15240/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/15240.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/16649/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/16649.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17072/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17072.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17160/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17160.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17167/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17167.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17175/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17175.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17369/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17369.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17671/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17671.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17922/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17922.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/17923/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/17923.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18055/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18055.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18121/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18121.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18134/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18134.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18153/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18153.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18241/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18241.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18266/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18266.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18282/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18282.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18355/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18355.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/15533_2/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/15533_2.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18432/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18432.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18559/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18559.gtf
#cp /c8000xd3/rnaseq-heath/Cufflinks/18694/transcripts.gtf /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/18694.gtf

cd /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare/
cuffcompare -p 8 -o ../Cuffmerge \
  -g /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
  -s /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa \
  15468.gtf 15533.gtf 15768.gtf 16286.gtf 16438.gtf 16488.gtf 16840.gtf 16929.gtf 16972.gtf \
  17049.gtf 17054.gtf 17068.gtf 17081.gtf 17087.gtf 17109.gtf 17115.gtf 17130.gtf 17701.gtf \
  17812.gtf 18294.gtf 18349.gtf 18596.gtf 18655.gtf 18687.gtf 15641.gtf 16024.gtf 16115.gtf \
  16385.gtf 16428.gtf 16491.gtf 16548.gtf 16810.gtf 16826.gtf 17048.gtf 17053.gtf 17071.gtf \
  17921-l1.gtf 15655.gtf 16483.gtf 16640.gtf 17013.gtf 17229.gtf 17333.gtf 17475.gtf \
  17543.gtf 17629.gtf 17835.gtf 18249.gtf 18372.gtf 18666.gtf 18983.gtf 19043.gtf 19052.gtf \
  A138.gtf A226.gtf 15240.gtf 16649.gtf 17072.gtf 17160.gtf 17167.gtf 17175.gtf 17369.gtf \
  17671.gtf 17922.gtf 17923.gtf 18055.gtf 18121.gtf 18134.gtf 18153.gtf 18241.gtf 18266.gtf \
  18282.gtf 18355.gtf 15533_2.gtf 18432.gtf 18559.gtf 18694.gtf 
echo "Finished Cuffmerge"
