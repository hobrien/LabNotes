- To download the sequences corresponding to a list of accessions:

    ```
    for seq in Accession Number list 
    do
       esearch -query $seq -db nuccore | efetch -format fasta
    done
    ```

