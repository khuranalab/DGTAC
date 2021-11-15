import pandas as pd
import re, sys
import os
import numpy as np
import random
from os import rename, listdir
import string


if len(sys.argv) != 3:
    print "Usage: python exp_splitCol.py [input adjusted tpm matrix] [count of common genes]"
    sys.exit(1)

tempPath=os.environ['PWD']
input = open(sys.argv[1], 'r')
commonGeneCnt = int(sys.argv[2])

exp_tr = pd.read_csv(input, index_col=0).transpose()

nSplit = 32

eachSplitGeneCnt = int(len(exp_tr.columns)/nSplit)+1
print eachSplitGeneCnt
for i in range(0, nSplit):
    if i < nSplit-1:
        splitExp = exp_tr.iloc[:,range(i*eachSplitGeneCnt, (i+1)*eachSplitGeneCnt)]
#        print i*eachSplitGeneCnt,(i+1)*eachSplitGeneCnt
        splitExp.to_csv(tempPath+'/.tmp/expAdjusted_split/expAdjusted_split_split'+str(i+1)+'.csv')
    else:
        splitExp = exp_tr.iloc[:,i*eachSplitGeneCnt:]
#        print i*eachSplitGeneCnt
        splitExp.to_csv(tempPath+'/.tmp/expAdjusted_split/expAdjusted_split_split'+str(i+1)+'.csv')

