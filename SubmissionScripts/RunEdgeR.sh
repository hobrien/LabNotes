# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 11 --max 20
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 11 --max 20 -a
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 12 --max 13
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 12 --max 13 -e 17046
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 13 --max 14
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 14 --max 15
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 15 --max 17
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 17 --max 20

# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 12 --max 15 --brainbank HDBRexpression
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 12 --max 13 --brainbank HDBRexpression
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 13 --max 14 --brainbank HDBRexpression
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 13 --max 15 --brainbank HDBRexpression
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 14 --max 15 --brainbank HDBRexpression
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 12 --max 15 --brainbank HDBRexpression -e 11373
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 12 --max 15 --brainbank HDBRexpression -e 11373,11775,11572
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 12 --max 13 --brainbank HDBRexpression -e 11373
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 13 --max 14 --brainbank HDBRexpression -e 11775,11572
# Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 13 --max 15 --brainbank HDBRexpression -e 11775,11572

for SampleID in 17812	18282	18349	18355	18687	19043	17921	19031	17221	18015	16428	16488	16640	16929	17053	17071	18529	18856
do
    Rscript ~/BTSync/FetalRNAseq/LabNotes/R/FvsMedgeR.r --min 14 --max 15 -e $SampleID
done

cp ~/BTSync/FetalRNAseq/Counts/MvsF_12_HDBRexpression_FDR_0.1_edgeR/tables/MaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/MaleUp12rep.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_12_HDBRexpression_FDR_0.1_edgeR/tables/FemaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/FemaleUp12rep.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_13_HDBRexpression_FDR_0.1_edgeR/tables/MaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/MaleUp13rep.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_13_HDBRexpression_FDR_0.1_edgeR/tables/FemaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/FemaleUp13rep.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_14_HDBRexpression_FDR_0.1_edgeR/tables/MaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/MaleUp14rep.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_14_HDBRexpression_FDR_0.1_edgeR/tables/FemaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/FemaleUp14rep.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_12_15_HDBRexpression_FDR_0.1_edgeR/tables/MaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/MaleUp12_14rep.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_12_15_HDBRexpression_FDR_0.1_edgeR/tables/FemaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/FemaleUp12_14rep.txt

cp ~/BTSync/FetalRNAseq/Counts/MvsF_13_15_HDBRexpression_FDR_0.1_edgeR/tables/FemaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/FemaleUp13_14rep.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_13_15_HDBRexpression_FDR_0.1_edgeR/tables/MaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/MaleUp13_14rep.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_13_15_HDBRexpression_FDR_0.1_edgeR/tables/MalevsFemale.background.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/BG13_14rep.txt


cp ~/BTSync/FetalRNAseq/Counts/MvsF_17_20_HDBR_excl_18432_FDR_0.1_edgeR/tables/FemaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/FemaleUp17_19.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_17_20_HDBR_excl_18432_FDR_0.1_edgeR/tables/MaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/MaleUp17_19.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_17_20_HDBR_excl_18432_FDR_0.1_edgeR/tables/MalevsFemale.background.txt ~/BTSync/FetalRNAseq/LabNotes/Results/BG17_19.txt

cp ~/BTSync/FetalRNAseq/Counts/MvsF_15_17_HDBR_excl_15641_FDR_0.1_edgeR/tables/FemaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/FemaleUp15_16.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_15_17_HDBR_excl_15641_FDR_0.1_edgeR/tables/MaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/MaleUp15_16.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_15_17_HDBR_excl_15641_FDR_0.1_edgeR/tables/MalevsFemale.background.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/BG15_16.txt

cp ~/BTSync/FetalRNAseq/Counts/MvsF_14_HDBR_FDR_0.1_edgeR/tables/FemaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/FemaleUp14.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_14_HDBR_FDR_0.1_edgeR/tables/MaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/MaleUp14.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_14_HDBR_FDR_0.1_edgeR/tables/MalevsFemale.background.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/BG14.txt

cp ~/BTSync/FetalRNAseq/Counts/MvsF_13_HDBR_excl_16491_FDR_0.1_edgeR/tables/FemaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/FemaleUp13.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_13_HDBR_excl_16491_FDR_0.1_edgeR/tables/MaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/MaleUp13.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_13_HDBR_excl_16491_FDR_0.1_edgeR/tables/MalevsFemale.background.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/BG13.txt

cp ~/BTSync/FetalRNAseq/Counts/MvsF_12_HDBR_FDR_0.1_edgeR/tables/FemaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/FemaleUp12.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_12_HDBR_FDR_0.1_edgeR/tables/MaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/MaleUp12.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_12_HDBR_FDR_0.1_edgeR/tables/MalevsFemale.background.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/BG12.txt

cp ~/BTSync/FetalRNAseq/Counts/MvsF_12_20_HDBR_excl_15641_18432_16491_PCW_FDR_0.1_edgeR/tables/FemaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/FemaleUp12_19.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_12_20_HDBR_excl_15641_18432_16491_PCW_FDR_0.1_edgeR/tables/MaleUp.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/MaleUp12_19.txt
cp ~/BTSync/FetalRNAseq/Counts/MvsF_12_20_HDBR_excl_15641_18432_16491_PCW_FDR_0.1_edgeR/tables/MalevsFemale.background.txt  ~/BTSync/FetalRNAseq/LabNotes/Results/BG12_19.txt

