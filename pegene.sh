#! /bin/sh

#################################################################
#################################################################
###       Protocol for pegene pipeline for additional samples ###
###                  -- Erica is a genius                     ###
#################################################################
#################################################################


#################################################################
# Input:
# - input sample list with cancer type as a second column (if info not available, please enter "Others" as cancer type)
# - Optional: sample source as the third column incase samples are from different batch/studies ('source' as column name).

sampleList=$1

# - peak log cpm matrix file, index=, column=
peakCpm=$2

# - bigwig file folders for each sample
bigwigFolder=$3

# - expression tpm matrix, index=ensembl gene id "ENSGXXXX", column=sample id
expTpm=$4

# - copy number matrix (optional, if non enter "none")
cnaMatrix=$5

# - panCancer model
panCancerModel=$6

mkdir $PWD/.tmp
mkdir $PWD/output

#################################################################
## Usage ########################################################
## Usage: pegene.sh sampleList.csv peakCpm.csv bigwigFolder expressionTpm.csv cnvMatrix.csv finalModel.sav
## Output prediction matrix will be output into your working directory
#################################################################

#################################################################
# Step 1. Combine expression files with TCGA 400 sample expression and batch correct
## - generate new batchInfo.txt file
## - formatting input expression tpm to log2(tpm+1)+1
## - batch correct merged expression, using TCGA 400 samples as ref.
## Output: /.tmp/expTpmMerged.csv; /.tmp/expTpmAdjusted.csv

python $PWD/scripts/expTpm2mergedExp.py $PWD/$expTpm $PWD/$sampleList
Rscript --vanilla $PWD/scripts/batchCorrect.R $PWD/.tmp/expTpmMerged.csv $PWD/.tmp/expTpmAdjusted.csv
commonGeneCnt=$(($(wc -l $PWD/.tmp/expTpmAdjusted.csv | cut -d' ' -f 1)-1))
mkdir $PWD/.tmp/expAdjusted_split
python $PWD/scripts/exp_splitCol.py $PWD/.tmp/expTpmAdjusted.csv $commonGeneCnt

#################################################################
# Step 2. Merge peak cpmMatrix.csv with ATAC-seq peak and batch correct

python $PWD/scripts/peakCpm2mergedPeak.py $PWD/$peakCpm
Rscript --vanilla $PWD/scripts/batchCorrect.R $PWD/.tmp/peakCpmMerged.csv $PWD/.tmp/peakCpmAdjusted.csv

#################################################################
# Step 3. Process the bigwig files to 5*100 bp windows, and process the same as in Step 2.
## - rename your bigwig files using sample names provided by sample list, and ended with ".bigwig"
## - bigwigFolder used in the initial version is "bigwigInsertCnt"

python $PWD/scripts/bigwig2bin.py $PWD/$bigwigFolder $PWD/$sampleList
for num in $(seq 1 5)
do
(
    python $PWD/scripts/binRdCnt2mergedCnt.py "$PWD/.tmp/bin${num}.csv" "${num}"
    Rscript --vanilla $PWD/scripts/batchCorrect_binNorm.R "$PWD/.tmp/binRdCntMerged_bin${num}.csv" "$PWD/.tmp/binRdCntAdjusted_bin${num}.csv"
) &
    sleep 5m
done
wait

#################################################################
# Step 4. (OPTIONAL) Calculate coefficient for all peaks around 1 Mb of a given region
## - To run faster, genes are splited into small files, in total 32 splits.
## - in total 6*32 jobs, peak, bin1~5
## - regress new data with source data in case the new data set has too little smaples

###### need to revise later ############## start ############
#for splitNum in $(seq 1 32)
#do
#    echo $splitNum
#    cp "calFeatureCoef_peakLog_splitExp.sh" "calFeatureCoef_peakLog_split${splitNum}.sh"
#    sed -i -e "s/Z/${splitNum}/g" "calFeatureCoef_peakLog_split${splitNum}.sh"
#
#    chmod +x "calFeatureCoef_peakLog_split${splitNum}.sh"
#    sbatch "calFeatureCoef_peakLog_split${splitNum}.sh"
#
#done
#
### - Merge splited coefficients in coefOutSplit into peakLog, bin1~5.
#
#./combineCoefOut_peakLog.sh
#./combineCoefOut_bin.sh
#python mergeCoef.py
###### need to revise later ############## end ############

#################################################################
# Step 5. Merge cnv file with ATAC-seq cnv and batch correct

python $PWD/scripts/cna2mergedCNA.py $PWD/$cnaMatrix
Rscript --vanilla $PWD/scripts/batchCorrect.R $PWD/.tmp/cnaMerged.csv $PWD/.tmp/cnaAdjusted.csv

wait
#################################################################
# Step 5. Calculate error term for each feature in each sample and the apply model
## - coefficient use either original coef from only TCGA samples, or use the newly calculated one from Step 4.
## - Also apply model
## - Output sample_DF for and mode application.
## - outputs are generated in output/df and output/predict


mkdir $PWD/output/df
mkdir $PWD/output/predict

sample_list=($(grep -vE '^sample' $PWD/$sampleList | cut -d',' -f 1))
echo $sample_list

for sample in "${sample_list[@]}"
do
(
    python $PWD/scripts/calError.py "${sample}"
    wait
    python $PWD/scripts/patientNetwork.py "${sample}"
)&
    sleep 10m
done

#################################################################
# Step 6. Using 0.5 cutoff to call peak-gene links in each patient
## - Output peak-gene patient matrix
## - Output cumulated count table

python $PWD/scripts/collectPeak2GeneInd.py 0.5



#################################################################
#################################################################
###       Following content are for downsteam analysis        ###
###           -- Not directly included in the pipeline        ###
#################################################################
#################################################################

#################################################################
# Step 7. Extract sample dataframe for peak-gene link of interests.
## - Output into linkSelect_df folder

mkdir linkSelect_df
sampleList=sampleList.txt
sample_list=($(grep -vE '^sample' $PWD/$sampleList | cut -d',' -f 1))
echo $sample_list

for sample in "${sample_list[@]}"
do
    (
    python $PWD/scripts/patientSelectLink_df.py "${sample}"
    )&
    sleep 10s
done
