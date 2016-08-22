find /c8000xd2/HiSeq4000-shared/heath/160817_K00267_0020_AHCWLGBBXX -name *.fastq.gz | grep -v Undetermined | xargs -n 1 qsub ../LabNotes/SubmissionScripts/FastQC.sh
find /c8000xd2/HiSeq4000-shared/heath/160817_K00267_0020_AHCWLGBBXX -name *.fastq.gz | grep -v Undetermined | sort | xargs -n 4 python ../LabNotes/Python/JoinReads.py | xargs -n 2 qsub ../LabNotes/SubmissionScripts/Tophat2.sh
