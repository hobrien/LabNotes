## Conservation of chloroplast development genes between Begonia and Dinoflagellates

The thylakoid stacking in the ocelloids of dinoflagellates is similar to the iridoplasts of *Begoina* {Gavelis:2015ir}. This has lead to the hypothesis that the same genetic programme is operating in both. Transcriptome data were generated in the ocelloid paper, but they are not publicly available (except for a couple of chloroplast genes), and are presumably of very low quality anyway. However, baring the possibility of HGT, the genes for thylakoid stacking would have to be conserved since the common ancestor of dinoflagellates and *Begoinia*. The concept of a "common ancestor" is complicated because there have been several rounds of ensymbiosis in the history of dinoflagellates. The genes that we are interested in have presumably been transfered from the nucleus of the symbiont to the nucleus of the host, but dinoflagellates have a complex mosaic of chloroplast development genes that have been transferred from green algae and red algae during different times in their history {Wisecaver:2011gf}. Nevertheless, it is probably worth looking for conservation between *Begoinia* and other species of dinoflagellate that do have whole-genome resources that are available. Unfortunately, dinoflagellate genomics seems to be a poorly developed field. A complete genome was published for the coral symbiont *Symbiodinium kawagutii* last year {Lin:2015ht} (along with a partial genome of *Symbiodinium minutum* from a few years back {Shoguchi:2013bx}). I've also been able to find a decent transcriptome assembly for *Lingulodinium polyedrum* {Beauchemin:2012hj}. There's supposedly also a transcriptome for *Alexandrium catenella* {Zhang:2014jf}, but I can only find the raw reads and not the assembly.

- I've downloaded protein sequences for *Symbiodinium kawagutii* (`~/Bioinformatics/Begonia/Dinos/SymbiodiniumAA.fa`) and transcriptome sequences for *Lingulodinium polyedrum* (`~/Bioinformatics/Begonia/Dinos/Lingulodinium.fa`)

- Make blast database for *Symbiodinium kawagutii*:
    - ```makeblastdb -in SymbiodiniumAA.fa -dbtype prot -out SymbiodiniumAA```
    - ```mkdir SymbiodiniumBlast```
- Write shell script to run BlastP all candidates:

```
    for file in `ls ../Candidates | grep Ath`
    do
        gene=$(echo $file | cut -d_ -f1)
        blastp -query ../Candidates/$file -db SymbiodiniumAA -evalue 1e-5 \
        -outfmt '6 qseqid qlen sacc slen pident length mismatch gapopen qstart qend qframe sstart send sframe evalue bitscore' \
        -out SymbiodiniumBlast/${gene}_Symbiodinium.bl
    done
```

- Make blast database for *Lingulodinium polyedrum*
    - ```makeblastdb -in Lingulodinium.fa -dbtype nucl -out Lingulodinium```
    - ```mkdir LingulodiniumBlast```
- Write shell script to run TblastN on all candidates:
```
    for file in `ls ../Candidates | grep Ath`
    do
        gene=$(echo $file | cut -d_ -f1)
        tblastn -query ../Candidates/$file -db Lingulodinium -evalue 1e-5 \ 
        -outfmt '6 qseqid qlen sacc slen pident length mismatch gapopen qstart qend qframe sstart send sframe evalue bitscore' \
        -out LingulodiniumBlast/${gene}_Lingulodinium.bl
    done 
```
- Remove empty blast results:

```find SymbiodiniumBlast -size  0 -print0 |xargs -0 -rm```
```find LingulodiniumBlast -size  0 -print0 |xargs -0 -rm```
         
- Create folders for sequences and alignments
    - ```mkdir SymbiodiniumSeqs```
    - ```mkdir SymbiodiniumAln```
    - ```mkdir LingulodiniumSeqs```
    - ```mkdir LingulodiniumAln```
    
- Write shell script to parse BlastP results:
```
    for file in `ls SymbiodiniumBlast`
    do
        gene=$(echo $file | cut -d_ -f1)
        ParseBlast.py -p blastp \
        --outfmt '6 qseqid qlen sacc slen pident length mismatch gapopen qstart qend qframe sstart send sframe evalue bitscore' \
        SymbiodiniumBlast/${gene}_Symbiodinium.bl SymbiodiniumAA.fa \
        > SymbiodiniumSeqs/${gene}_SymbiodiniumAA.fa
        mafft --add SymbiodiniumSeqs/${gene}_SymbiodiniumAA.fa \
        ../Candidates/${gene}_aln.fa > SymbiodiniumAln/${gene}_SymbiodiniumAA_aln.fa 
    done 
```

- Notes:
    - Had to exclude a few APG6 genes to get a comprehensible alignment. The removed sequences are in `~/Bioinformatics/Begonia/Candidates/APG6_excluded_seqs.fa`
    - The alignment of the AS1 homolog is quite crap. The Evalue of this hit was 2e-07
    - FZL is also pretty bad. This one is 7e-08
    - THF1 alignment isn't great. Evalue is 9e-17, so it probably is legit
    
- Write shell script to parse TblastN results:
```
    for file in `ls LingulodiniumBlast`
    do
        gene=$(echo $file | cut -d_ -f1)
        ParseBlast.py -p tblastn -t \
        --outfmt '6 qseqid qlen sacc slen pident length mismatch gapopen qstart qend qframe sstart send sframe evalue bitscore' \
        LingulodiniumBlast/${gene}_Lingulodinium.bl Lingulodinium.fa \
        > LingulodiniumSeqs/${gene}_Lingulodinium.fa
        mafft --add LingulodiniumSeqs/${gene}_Lingulodinium.fa \s
        ../Candidates/${gene}_aln.fa > LingulodiniumAln/${gene}_LingulodiniumAA_aln.fa 
    done 
```

- Notes:
    - Had to exclude a few ARC3 genes to get a comprehensible alignment. The removed sequences are in `~/Bioinformatics/Begonia/Candidates/ARC3_excluded_seqs.fa`
   - The alignment of the AS1 homolog is quite crap. The Evalue of this hit was 3e-07
   - The ATTERC homolog is just the TerC superfamily domain (Evalue is 3e-48)
   - The COP1 homolog is just the WD40 superfamily domain. The sequence seems very diverged from the plant homologs and the sequence is most similar to peroxisomal biogenesis factor 7 (Evalue is 4e-18)
   - The alignment of the PYG7 homolog is quite crap. The Evalue of this hit was 2e-08

- Write shell scripts to extract Evalue and name of top hit for each gene 
```

    for file in `ls SymbiodiniumBlast`
    do
        gene=$(echo $file | cut -d_ -f1)
        stats=$(head -1 SymbiodiniumBlast/$file | cut -f3,15)
        echo $gene $stats 
    done 
```

``` bash shell_script.sh > Symbiodinium.txt```

```

    for file in `ls LingulodiniumBlast`
    do
        gene=$(echo $file | cut -d_ -f1)
        stats=$(head -1 LingulodiniumBlast/$file | cut -f3,15)
        echo $gene $stats 
    done 
```

``` bash shell_script.sh > Lingulodinium.txt```

- Write shell scripts to extract the sequence of the top hit from each genome, blast it against swissprot, and write the subject title  (only available in blast+ version 2.2.27 or higher) to file

```mkdir SymbiodiniumSwissProt```

```

    for file in `ls SymbiodiniumBlast`
    do
        gene=$(echo $file | cut -d_ -f1)
        head -1 SymbiodiniumBlast/$file > temp.bl
        ParseBlast.py -p blastp  \
        --outfmt '6 qseqid qlen sacc slen pident length mismatch gapopen qstart qend qframe sstart send sframe evalue bitscore' \
        temp.bl SymbiodiniumAA.fa > temp.fa
        blastp -remote -query temp.fa -db swissprot -max_target_seqs 1 \
        -outfmt '6 qseqid sacc stitle evalue' | head -1 \
        > SymbiodiniumSwissProt/${gene}_Symbiodinium.bl
    done 
    
```

```mkdir LingulodiniumSwissProt```

```

    for file in `ls LingulodiniumBlast`
    do
        gene=$(echo $file | cut -d_ -f1)
        head -1 LingulodiniumBlast/$file > temp.bl
        ParseBlast.py -p tblastn -t  \
        --outfmt '6 qseqid qlen sacc slen pident length mismatch gapopen qstart qend qframe sstart send sframe evalue bitscore' \
        temp.bl Lingulodinium.fa > temp.fa
        blastp -remote -query temp.fa -db swissprot -max_target_seqs 1 \
        -outfmt '6 qseqid sacc stitle evalue' | head -1 \
        > LingulodiniumSwissProt/${gene}_Lingulodinium.bl
    done 
    
```

```    rm temp.bl temp.fa ```

- I can only do about 6 or 8 remote blast searches before they stop working. If I use a VPN, they start working again, so this is clearly an NCBI-imposed limit. I can't find any documentation about this, but it is pretty annoying.

- Wrote a [doit](http://pydoit.org) script to handle blasting without having to redo completed steps every time the task is rerun:
    - ```doit -n 3 -f BlastCandidates.py```
    - currently, this only works with the Lingulodinium genes. The Symbiodinium ones are all done anyway so I'm not sure it's worth redoing them as part of this script (jobs that aren't done by doit have to be repeated because they're not in the doit database)
    
- The rest of this is not going to go in the doit script because of the problem of redoing all the blast searches.
- Concatinate top hits and merge with original blast results:
    - ```cat SymbiodiniumSwissProt/* >Symbiodinium_top_hits.txt```
- Use dplyr in R to merge:
    ```
    library(dplyr)
    
Symbiodinium <- read.table("~/Bioinformatics/Begonia/Dinos/Symbiodinium.txt", quote="\"", comment.char="")
Symbiodinium_top_hits <- read.delim("~/Bioinformatics/Begonia/Dinos/Symbiodinium_top_hits.txt", header=FALSE)
Symbiodinium <- full_join(Symbiodinium, Symbiodinium_top_hits, by=c("V2" = "V1"))
names(Symbiodinium)<- c('gene', 'orf', 'Evalue', 'top_hit', 'top_hit_description', 'top_hit_evalue')
write.table(Symbiodinium, file="~/Bioinformatics/Begonia/Dinos/Symbiodinium_merge.txt", sep="\t", quote=F, row.names=F)

Lingulodinium <- read.table("~/Bioinformatics/Begonia/Dinos/Lingulodinium.txt", quote="\"", comment.char="")
Lingulodinium_top_hits <- read.delim("~/Bioinformatics/Begonia/Dinos/Lingulodinium_top_hits.txt", header=FALSE)
Lingulodinium <- full_join(Lingulodinium, Lingulodinium_top_hits, by=c("V2" = "V1"))
names(Lingulodinium)<- c('gene', 'orf', 'Evalue', 'top_hit', 'top_hit_description', 'top_hit_evalue')
write.table(Lingulodinium, file="~/Bioinformatics/Begonia/Dinos/Lingulodinium_merge.txt", sep="\t", quote=F, row.names=F)
```

- This is looking good, but it would be nice to know the organism for each top hit
    - This will extract this info for an accession number:
        - ```efetch -db protein -id Q8S0J7 -format gpc |xtract -insd INSDSeq_organism```
        
    - This will do it for all:
        ```
        for acc in `cut -f2 Symbiodinium_top_hits.txt`
     do
            efetch -db protein -id $acc -format gpc |xtract -insd INSDSeq_organism >> Symbiodinium_top_hit_taxon.txt
    done
    for acc in `cut -f2 Lingulodinium_top_hits.txt`
    do
            efetch -db protein -id $acc -format gpc |xtract -insd INSDSeq_organism >> Lingulodinium_top_hit_taxon.txt
    done
```    
