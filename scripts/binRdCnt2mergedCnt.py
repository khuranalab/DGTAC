
import pandas as pd
import re, sys
import os
import numpy as np
import random
import scipy.stats as stats
from scipy.stats import norm


if len(sys.argv) != 3:
    print "Usage: python binRdCnt2mergedCnt [input bin matrix] [bin number]"
    sys.exit(1)

tempPath=os.environ['PWD']
input = open(sys.argv[1], 'r')
binN = str(sys.argv[2])
#sampleInfo = pd.read_csv(open(sys.argv[2], 'r'), index_col=0)
tcgaSampleInfo = pd.read_csv(tempPath+'/source/tcgaSampleInfo.csv', index_col=0)

rdCnt = pd.read_csv(input, index_col=0).transpose()
#print rdCnt.head(5)
#print rdCnt.isnull().any().any()

# TCGA sample peak matrix, index=sample ID, columns=gene names
tcgaCpmLog = pd.read_csv(tempPath+'/source/bin'+binN+'_tsp.csv', index_col=0)
#print tcgaCpmLog.isnull().any().any()

# check two peak files are consistent
if len(set(tcgaCpmLog.index)&set(rdCnt.index)) != len(tcgaCpmLog.index):
    print 'Error: two peak sets don\'t match, please regenerate your count Matrix'

else:
    # To keep the batchInfo consistent with expression batch info
    batchInfo = pd.read_csv(tempPath+'/.tmp/batchInfo.csv', index_col=0)

    # Output merged expression dataset.
    merged4combat = pd.concat([tcgaCpmLog, rdCnt], axis=1, sort=False)
    merged4combat = merged4combat[batchInfo.index.tolist()]

    if ~merged4combat.isnull().any().any():
        merged4combat.to_csv(tempPath+'/.tmp/binRdCntMerged_bin'+binN+'.csv')
    else:
        print 'Error: null values in the merged bin read count table.'
