
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
import time
import datetime

import random


## temp path on node
tempPath = os.environ['PWD']
selectedLink = pd.read_csv(tempPath+'/ESR1_selectedLinkList.txt', index_col=0, header=None).index.tolist()

eachPatient = str(sys.argv[1])
patientDF_test0 = pd.read_csv(tempPath+'/output/df/'+eachPatient+'_df.csv', index_col=0)

print 'All pairs', len(patientDF_test0.index)

#patientDF_test0 = patientDF_test0.fillna(0)
#patientDF_test = patientDF_test0
#for eachError in patientDF_test0.columns[9:]:
#    #     print patientDF_test[(abs(patientDF_test[eachError]) <= 1)]
#    patientDF_test = patientDF_test[(abs(patientDF_test[eachError]) < 10)]
#
#print 'Pairs with error < 10, aka gene expression = 1', len(patientDF_test)
#patientDF_test.head()

patientDF_test0.loc[selectedLink,:].to_csv(tempPath+'/linkSelect_df/selectLink_'+eachPatient+'.csv')
