
# coding: utf-8

# # Gather information for each patient for predicting network
# ### peak/bin, peak error/bin error

# In[1]:


import sys
import re
import pandas as pd
import numpy as np
import os
import time
import datetime


# ## Merge coefOut into coefOut_merged

# ### Features used

# In[3]:
tempPath = os.environ['PWD']

cutoff = float(sys.argv[1])

batchInfo = pd.read_csv(tempPath+'/.tmp/batchInfo.csv', index_col=0)
patientList = batchInfo[batchInfo['batch'] != 1].index.tolist()


n=0
peak2GeneCntDic = {}
cumulatedCnt = []

for eachPatient in patientList:
    print eachPatient
    yPredict_tmp = pd.read_csv(tempPath+'/output/predict/yPredict_'+eachPatient+'.csv', header=None, index_col=0)

    print 'all pairs:', len(yPredict_tmp.index)
    yPredict_tmp.columns = ['prob','label']
    yPredict_tmp = yPredict_tmp[yPredict_tmp['prob'] > cutoff]
    print 'true pairs:', len(yPredict_tmp.index)

    for eachPair in yPredict_tmp.index:
        peak2GeneCntDic.setdefault(eachPair, [0]*len(patientList))
        peak2GeneCntDic[eachPair][n] += 1

#    if n == 0:
#        peak2GeneDF['cnt'] = [0]*len(yPredict_tmp.index)
#        peak2GeneDF.index = yPredict_tmp.index


#    for eachPair in yPredict_tmp.index:
#        try:
#            peak2GeneDF.at[eachPair, 'cnt'] += 1
#        except:
#            addIn = pd.DataFrame([[1, np.nan]], columns=['cnt', 'freq'], index=[eachPair])
#            peak2GeneDF = pd.concat([peak2GeneDF,addIn])


#    print len(peak2GeneDF.index)

    print 'total true pairs:',len(peak2GeneCntDic.keys())
    cumulatedCnt.append(len(peak2GeneCntDic.keys()))
    n+=1

peak2GeneMatrix = pd.DataFrame.from_dict(peak2GeneCntDic, orient='index', columns=patientList)
peak2GeneMatrix.to_csv(tempPath+'/output/peak2GeneIndMatrix_'+str(cutoff)+'.csv')
#peak2GeneDF.to_csv('peak2GeneCnt.csv')


cumCntOut = open(tempPath+'/output/cumulatedCnt_'+str(cutoff)+'.txt','w')

for cnt in cumulatedCnt:
    cumCntOut.write(str(cnt)+'\n')

cumCntOut.close()
