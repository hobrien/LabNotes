#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:$PATH

bash ~/LabNotes/SubmissionScripts/mpileup.sh chr2:184936178 #ZNF804A
bash ~/LabNotes/SubmissionScripts/mpileup.sh chr1:98046571 #MIR137HG
bash ~/LabNotes/SubmissionScripts/mpileup.sh chr2:192775750 #PCGEM1
bash ~/LabNotes/SubmissionScripts/mpileup.sh chr6:30186422 #TRIM26
bash ~/LabNotes/SubmissionScripts/mpileup.sh chr8:2937908 #CSMD1
bash ~/LabNotes/SubmissionScripts/mpileup.sh chr8:88039069 #MMP16
bash ~/LabNotes/SubmissionScripts/mpileup.sh chr10:103054405 #CNNM2
bash ~/LabNotes/SubmissionScripts/mpileup.sh chr10:103089711 #NT5C2
bash ~/LabNotes/SubmissionScripts/mpileup.sh chr11:125609468 #STT3A
bash ~/LabNotes/SubmissionScripts/mpileup.sh chr18:54903114 #CCDC68
bash ~/LabNotes/SubmissionScripts/mpileup.sh chr18:55228300 #TCF4

