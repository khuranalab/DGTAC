
# coding: utf-8

# # Read the TF-peak matrix, and build the network for patient

# In[1]:


#Importing libraries.

import re
import sys
import numpy as np
import pandas as pd
import os
from scipy import interp
import scipy.stats as stats
import pickle
import time
import datetime

import random
from math import log

import networkx as nx

from networkx.algorithms import bipartite

print nx.__version__


## temp path on node
tempPath = os.environ['PWD']
if len(sys.argv) != 2:
    print "Usage: python patientNetwork [sample]"
    sys.exit(1)

rf = pickle.load(open(tempPath+'/source/finalModel_mcf7_balancedSet_all.sav', 'rb'))

eachPatient = str(sys.argv[1])
patientDF_test0 = pd.read_csv(tempPath+'/output/df/'+eachPatient+'_df.csv', index_col=0)

print 'All pairs', len(patientDF_test0.index)

patientDF_test0 = patientDF_test0.fillna(0)
patientDF_test = patientDF_test0
for eachError in patientDF_test0.columns[9:]:
    #     print patientDF_test[(abs(patientDF_test[eachError]) <= 1)]
    patientDF_test = patientDF_test[(abs(patientDF_test[eachError]) < 10)]

print 'Pairs with error < 10, aka gene expression = 1', len(patientDF_test)
patientDF_test.head()
# In[3]:

yPredict = rf.predict(patientDF_test)
yProb = (rf.predict_proba(patientDF_test)[:,1])

yPredictDic = {}
yPredictOutput = open(tempPath+'/output/predict/yPredict_'+eachPatient+'.csv', 'w')
yProbList = list(yProb)
for n in range(0,len(patientDF_test.index)):
    yPredictDic.setdefault(patientDF_test.index[n], yProbList[n])
    yPredictDic[patientDF_test.index[n]] = yProbList[n]
    yPredictOutput.write(patientDF_test.index[n]+','+str(yProbList[n])+','+yPredict[n]+'\n')

print len(yPredictDic.keys())
yPredictOutput.close()



# In[10]:


#outputGraph = open('patientNetwork/TCGA-EJ-7318.netwotk', 'w')
#nx.write_gml(egNetwrkG, outputGraph)
#
#outputGraph.close()

