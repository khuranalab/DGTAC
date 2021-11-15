
import pandas as pd
import re, sys
import os
import numpy as np
import random
import scipy.stats as stats
from scipy.stats import norm


if len(sys.argv) != 3:
    print "Usage: python expTpm2mergedExp.py [input tpm matrix] [sample info table]"
    sys.exit(1)

tempPath=os.environ['PWD']
input = open(sys.argv[1], 'r')
sampleInfo = pd.read_csv(open(sys.argv[2], 'r'), index_col=0)
tcgaSampleInfo = pd.read_csv(tempPath+'/source/tcgaSampleInfo.csv', index_col=0)

ensemblTb = pd.read_csv(tempPath+'/source/Homo_sapiens.GRCh38.79.geneID', header=None, sep='\t', index_col=4)

tpm = pd.read_csv(input, index_col=0)
tpm = tpm.loc[~(tpm==0).all(axis=1)]

geneName = []
for geneID in tpm.index:
    try:
        geneName.append(ensemblTb.at[geneID,5])
    except:
        geneName.append(geneID)
        print(geneID)
        continue

# For gene names repetes, keep the first one.
tpm.index = geneName
tpm = tpm.loc[~tpm.index.duplicated(keep='first')]

# log2 transform tpm, as well as add artifitial 1 for later use (avoiding denominater 0).
tpm = tpm.apply(lambda x: np.log2(x+1)+1)

# TCGA sample expression matrix, index=sample ID, columns=gene names
tcgaTpmLog = pd.read_csv(tempPath+'/source/selectedExpression_log2tpm.csv', index_col=0)

# find common genes
commonGenes = list(set(tcgaTpmLog.columns)&set(tpm.index))
print 'Shared genes names between TCGA RNA-seq and the new dataset: ', len(commonGenes)
# Set environment variables commonGeneCnt

# Output merged expression dataset.
merged4combat = pd.concat([tcgaTpmLog.loc[:,commonGenes].transpose(),tpm.loc[commonGenes,:]], axis=1, sort=False)
# Output merged batch info dataset.

if 'source' in sampleInfo.columns:
    sourceList = list(set(sampleInfo['source']))
    sample_batch_list = []
    for eachSample in sampleInfo.index:
        sample_batch_list.append(sourceList.index(sampleInfo.loc[eachSample,'source'])+2)
        
    tcgaSampleInfo.index.name = sampleInfo.index.name
    merged4combatInfo = pd.concat([tcgaSampleInfo, sampleInfo.loc[:,['cancerType']]], axis=0, sort=False)
    print(merged4combatInfo.head(2))
    merged4combatInfo['batch'] = [1]*len(tcgaSampleInfo.index)+sample_batch_list
else:
    merged4combatInfo = pd.concat([tcgaSampleInfo, sampleInfo], axis=0, sort=False)
    merged4combatInfo['batch'] = [1]*len(tcgaSampleInfo.index)+[2]*len(sampleInfo.index)

merged4combatInfo.to_csv(tempPath+'/.tmp/batchInfo.csv')

if ~merged4combat.isnull().any().any():
    merged4combat.to_csv(tempPath+'/.tmp/expTpmMerged.csv')
else:
    print 'Error: null values in the merged expression table.'


