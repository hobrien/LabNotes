## Conservation of chloroplast development genes between Begonia and Dinoflagellates

The thylakoid stacking in the ocelloids of dinoflagellates is similar to the iridoplasts of *Begoina* {Gavelis:2015ir}. This has lead to the hypothesis that the same genetic programme is operating in both. Transcriptome data were generated in the ocelloid paper, but they are not publicly available (except for a couple of chloroplast genes), and are presumably of very low quality anyway. However, baring the possibility of HGT, the genes for thylakoid stacking would have to be conserved since the common ancestor of dinoflagellates and *Begoinia*. The concept of a "common ancestor" is complicated because there have been several rounds of ensymbiosis in the history of dinoflagellates. The genes that we are interested in have presumably been transfered from the nucleus of the symbiont to the nucleus of the host, but dinoflagellates have a complex mosaic of chloroplast development genes that have been transferred from green algae and red algae during different times in their history {Wisecaver:2011gf}. Nevertheless, it is probably worth looking for conservation between *Begoinia* and other species of dinoflagellate that do have whole-genome resources that are available. Unfortunately, dinoflagellate genomics seems to be a poorly developed field. A complete genome was published for the coral symbiont *Symbiodinium kawagutii* last year {Lin:2015ht} (along with a partial genome of *Symbiodinium minutum* from a few years back {Shoguchi:2013bx}). I've also been able to find a decent transcriptome assembly for *Lingulodinium polyedrum* {Beauchemin:2012hj}. There's supposedly also a transcriptome for *Alexandrium catenella* {Zhang:2014jf}, but I can only find the raw reads and not the assembly.

- I've downloaded protein sequences for *Symbiodinium kawagutii* (Bioinformatics/Begonia/Dinos/SymbiodiniumAA.fa) and transcriptome sequences for *Lingulodinium polyedrum* (Bioinformatics/Begonia/Dinos/Lingulodinium.fa).

- Make blast database for *Symbiodinium kawagutii*:
    - ```makeblastdb -in SymbiodiniumAA.fa -dbtype prot -out SymbiodiniumAA```
    - ```mkdir SymbiodiniumBlast```
- Write shell script to run BlastP all candidates:
    - ```for file in `ls ../Candidates | grep Ath`
             do
             gene=$(echo $file | cut -d_ -f1)
             blastp -query ../Candidates/$file -db SymbiodiniumAA -evalue 1e-5 -outfmt '6                           qseqid qlen sacc slen pident length mismatch gapopen qstart qend qframe sstart send sframe evalue bitscore' -out SymbiodiniumBlast/${gene}_Symbiodinium.bl
             done 
      ```
- Make blast database for *Lingulodinium polyedrum*
    - ```makeblastdb -in Lingulodinium.fa -dbtype nucl -out Lingulodinium```
    - ```mkdir LingulodiniumBlast```
- Write shell script to run TblastN on all candidates:
    - ```for file in `ls ../Candidates | grep Ath`
             do
             gene=$(echo $file | cut -d_ -f1)
             tblastn -query ../Candidates/$file -db Lingulodinium -evalue 1e-5 -outfmt '6                           qseqid qlen sacc slen pident length mismatch gapopen qstart qend qframe sstart send sframe evalue bitscore' -out LingulodiniumBlast/${gene}_Lingulodinium.bl
             done 
      ```
