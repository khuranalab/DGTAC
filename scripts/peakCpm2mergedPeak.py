
import pandas as pd
import re, sys
import os
import numpy as np
import random
import scipy.stats as stats
from scipy.stats import norm


if len(sys.argv) != 2:
    print "Usage: python peakCpm2mergedPeak.py [input cpm matrix]"
    sys.exit(1)

tempPath=os.environ['PWD']
input = open(sys.argv[1], 'r')
#sampleInfo = pd.read_csv(open(sys.argv[2], 'r'), index_col=0)
tcgaSampleInfo = pd.read_csv(tempPath+'/source/tcgaSampleInfo.csv', index_col=0)

cpm = pd.read_csv(input, index_col=0)

# TCGA sample peak matrix, index=sample ID, columns=gene names
tcgaCpmLog = pd.read_csv(tempPath+'/source/selectedPeakLog.csv', index_col=0)

# check two peak files are consistent
if len(tcgaCpmLog.index) != len(cpm.index):
    print 'Error: two peak sets don\'t match, please regenerate your cpm Matrix'

else:
    # To keep the batchInfo consistent with expression batch info
    batchInfo = pd.read_csv(tempPath+'/.tmp/batchInfo.csv', index_col=0)

    # Output merged expression dataset.
    merged4combat = pd.concat([tcgaCpmLog, cpm], axis=1, sort=False)
    merged4combat = merged4combat[batchInfo.index.tolist()]

    if ~merged4combat.isnull().any().any():
        merged4combat.to_csv(tempPath+'/.tmp/peakCpmMerged.csv')
    else:
        print 'Error: null values in the merged expression table.'
