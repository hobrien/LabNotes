- [This](https://sequencing.qcfail.com/articles/genomic-sequence-not-in-the-genome-assembly-creates-mapping-artefacts/) post points out that including as much of the reference as possible during mapping can improve analyses even if the extra sequence is filtered out downstream because it prevents the mapper from trying to map reads from these extra regions to the main chromosomes

- There are a bunch of different annotations. I need to figure out which I am going to use:
    - [GENCODE](http://www.gencodegenes.org/):
        - [hg38 BED](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_GENCODE_v23.bed.gz/download)
            - might be safer to [convert](http://bedops.readthedocs.org/en/latest/content/reference/file-management/conversion/gtf2bed.html) the GTF igenome GTF to BED to ensure that the formatting is identical
        - GTF: NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf
            - this is version version 23 (Ensembl 81). Current version is 24
            - this appears to correspond to GENCODE Comprehensive ALL, but the names are different
    - [UCSC](https://genome.ucsc.edu/):
        - [hg38 igenome](ftp://ftp.illumina.com/Homo_sapiens/UCSC/hg38)
        - [hg38 knownGene BED file](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_UCSC_knownGene.bed.gz/download)
            - the knowGene GTF doesn't work with BAMQC because it only has exon features
    - [NCBI](http://www.ncbi.nlm.nih.gov/refseq/):
        - [hg38 igenome (with decoy sequences)](ftp://ftp.illumina.com/Homo_sapiens/NCBI/GRCh38Decoy)
        - [RefSeq BED](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz/download)?
    - [Ensembl](http://www.ensembl.org/index.html):
        - [hg19 igenome](ftp://ftp.illumina.com/Homo_sapiens/Ensembl/GRCh37)
    - [rRNA sequences](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/GRCh38_rRNA.bed.gz/download)
    - "the terms 'Ensembl annotation' and 'GENCODE annotation' are thus synonymous when referring to human" {Frankish:2015iu}
    - "the GENCODE Comprehensive set is larger and represents more transcriptional complexity than RefSeq" {Frankish:2015iu}
    
- [Differences between GTF and GFF3](http://blog.nextgenetics.net/?e=27)

- Based on Google and [this](http://biorxiv.org/content/biorxiv/early/2013/11/22/000851.full.pdf) paper, it look slike STAR supports clipping of reads but Tophat doesn't. Might be an argument for quality trimming if I'm using tophat

- [Brainspan](http://www.brainspan.org/rnaseq/search/index.html): brain microarray data during development