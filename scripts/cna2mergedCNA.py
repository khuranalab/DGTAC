
import pandas as pd
import re, sys
import os
import numpy as np
import random
import scipy.stats as stats
from scipy.stats import norm


if len(sys.argv) != 2:
    print "Usage: python cna2mergedCNA.py [input cna matrix]"
    sys.exit(1)

tempPath=os.environ['PWD']
input = open(sys.argv[1], 'r')
#sampleInfo = pd.read_csv(open(sys.argv[2], 'r'), index_col=0)
tcgaSampleInfo = pd.read_csv(tempPath+'/source/tcgaSampleInfo.csv', index_col=0)

cpm = pd.read_csv(input, index_col=0)
#print cpm.isnull().any().any()
# TCGA sample peak matrix, index=sample ID, columns=gene names
tcgaCpmLog = pd.read_csv(tempPath+'/source/selectedCNA.csv', index_col=0)
#print tcgaCpmLog.isnull().any().any()

# get common genes:
commonGene=list(set(tcgaCpmLog.index)&set(cpm.index))
print len(tcgaCpmLog.index), len(cpm.index), len(commonGene)

# To keep the batchInfo consistent with expression batch info
batchInfo = pd.read_csv(tempPath+'/.tmp/batchInfo.csv', index_col=0)

# Output merged expression dataset.
merged4combat = pd.concat([tcgaCpmLog.loc[commonGene,:], cpm.loc[commonGene,:]], axis=1, sort=False)
merged4combat = merged4combat[batchInfo.index.tolist()]

merged4combat = merged4combat.fillna(0)
print len(merged4combat.index)

if ~merged4combat.isnull().any().any():
    merged4combat.to_csv(tempPath+'/.tmp/cnaMerged.csv')
else:
    print 'Error: null values in the merged expression table.'
