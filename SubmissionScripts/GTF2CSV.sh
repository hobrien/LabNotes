for input in `find /c8000xd3/rnaseq-heath/Cufflinks/Cuffcompare -name "*.gtf" | grep -v sorted`
do
    file=${input##*/}
    sampleID=${file%%.*} 
    cat $input |python ~/LabNotes/Python/GTF2CSV.py $sampleID >> /c8000xd3/rnaseq-heath/Cufflinks/AllGTF.csv 
done
