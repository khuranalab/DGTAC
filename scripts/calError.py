
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

## temp path on node
if len(sys.argv) != 2:
    print "Usage: python calError.py [sample]"
    sys.exit(1)

tempPath=os.environ['PWD']
eachPatient = str(sys.argv[1])

### Read in coef output
coefOut_mergedCol = ['peakLog_intercept','peakLog_coef','bin1_intercept','bin1_coef','bin2_intercept','bin2_coef','bin3_intercept','bin3_coef','bin4_intercept','bin4_coef','bin5_intercept','bin5_coef']
coefOut_merged_dtypeDic = {}
for pt in coefOut_mergedCol:
    coefOut_merged_dtypeDic.setdefault(pt, np.float64)
coefOut_merged = pd.read_csv(tempPath+'/source/coefOut_merged.csv', index_col=0, dtype=coefOut_merged_dtypeDic)

featureList = ['peakLog','bin1','bin2','bin3','bin4','bin5']

### Read in batch corrected expression
expression = pd.DataFrame(pd.read_csv(tempPath+'/.tmp/expTpmAdjusted.csv', index_col=0).loc[:,eachPatient])
print expression.head(2)

### Read in batch corrected peak log
peak = pd.DataFrame(pd.read_csv(tempPath+'/.tmp/peakCpmAdjusted.csv', index_col=0).loc[:,eachPatient])
print peak.head(2)

### Read in bins
bin1 = pd.DataFrame(pd.read_csv(tempPath+'/.tmp/binRdCntAdjusted_bin1.csv', index_col=0).loc[:, eachPatient])
bin2 = pd.DataFrame(pd.read_csv(tempPath+'/.tmp/binRdCntAdjusted_bin2.csv', index_col=0).loc[:, eachPatient])
bin3 = pd.DataFrame(pd.read_csv(tempPath+'/.tmp/binRdCntAdjusted_bin3.csv', index_col=0).loc[:, eachPatient])
bin4 = pd.DataFrame(pd.read_csv(tempPath+'/.tmp/binRdCntAdjusted_bin4.csv', index_col=0).loc[:, eachPatient])
bin5 = pd.DataFrame(pd.read_csv(tempPath+'/.tmp/binRdCntAdjusted_bin5.csv', index_col=0).loc[:, eachPatient])
print bin1.head(2)

####### Add CNA later
cna = pd.read_csv(tempPath+'/.tmp/cnaAdjusted.csv', index_col=0)
cna_tr = cna.transpose()

# ## Gene Location
input = open(tempPath+'/source/refSeq_geneNames_hg38.txt', 'r')
geneTssDic = {}
for line in input:
    line = line.strip()
    if line.startswith('chrom'):
        continue
    fields = line.split('\t')
    
    geneTssDic.setdefault(fields[3], [])
    if fields[4] == '+':
        geneTssDic[fields[3]].append((fields[0],int(fields[1])))
    else:
        geneTssDic[fields[3]].append((fields[0],int(fields[2])))


print len(geneTssDic)
print geneTssDic['MUC7']


# ## Peak Location

# In[23]:


pkLoc = pd.read_csv(tempPath+'/source/peakLocation.csv', index_col=0)

panCanPeakLoc = pkLoc.iloc[:,[0,1,2]]
print panCanPeakLoc.head(2)


# # Calculate errors for eachFeatures and each patient
print eachPatient, ' start'
startTime = time.time()
patientOut = open(tempPath+'/output/df/'+eachPatient+'_df.csv', 'w')
patientOut.write('gene-peak,expression,disTSS,cna,'+','.join(featureList)+','+','.join(map(lambda feature: feature+'_error', featureList))+'\n')

for eachPair in coefOut_merged.index:
    g = re.match('(.+)-(ACC|ESCA|GBM|PCPG|STAD|UCEC|THCA|CESC|LIHC|CHOL|HNSC|SKCM|COAD|BLCA|TGCT|LUSC|MESO|KIRP|LGG|PRAD|LUAD|BRCA|KIRC)_(.+)', eachPair)
    pGene = g.group(1)
    pPeak = g.group(2)+'_'+g.group(3)
    
    if not pGene in geneTssDic.keys():
        continue
    
    
    patientOut.write(eachPair)
    
    try:
        yExpr = float(expression.loc[pGene, eachPatient])
    except:
        yExpr = 1

    patientOut.write(','+str(yExpr))
    
    
    tssList = []
    for i in geneTssDic[pGene]:
        tssList.append(int(i[1]))

    #         print brcaPeakLoc.loc[pPeak][2], brcaPeakLoc.loc[pPeak][1]
    #         print tssList
    peakCenter = (panCanPeakLoc.loc[pPeak][2]+panCanPeakLoc.loc[pPeak][1])/2
    distance2tss = abs(peakCenter-np.mean(tssList))
    patientOut.write(','+str(distance2tss))

    if pGene in cna.index:
        geneCNA = cna.loc[pGene, eachPatient]
    else:
        geneCNA = 0
    patientOut.write(','+str(geneCNA))

    patientOut.write(','+str(peak.loc[pPeak, eachPatient]))
    patientOut.write(','+str(bin1.loc[pPeak, eachPatient]))
    patientOut.write(','+str(bin2.loc[pPeak, eachPatient]))
    patientOut.write(','+str(bin3.loc[pPeak, eachPatient]))
    patientOut.write(','+str(bin4.loc[pPeak, eachPatient]))
    patientOut.write(','+str(bin5.loc[pPeak, eachPatient]))
    
    featureValues = [peak.loc[pPeak, eachPatient],bin1.loc[pPeak, eachPatient],bin2.loc[pPeak, eachPatient],bin3.loc[pPeak, eachPatient],bin4.loc[pPeak, eachPatient],bin5.loc[pPeak, eachPatient]]

    for colN in range(0, len(featureList)):
        ###### Error calculated as in Cao et al.
        if yExpr == 1:
            error = 10
        elif coefOut_merged.ix[eachPair, 2*colN+1] == 0:
            error = 10
#            yPredicExpr = coefOut_merged.ix[eachPair, 2*colN] + coefOut_merged.ix[eachPair, 2*colN+1]*float(featureValues[colN])
#            error = (yPredicExpr - yExpr)/abs(yExpr)
#            print eachPair, yExpr, yPredicExpr, error
        else:
            yPredicExpr = coefOut_merged.ix[eachPair, 2*colN] + coefOut_merged.ix[eachPair, 2*colN+1]*float(featureValues[colN])
            error = (yPredicExpr - yExpr)/abs(yExpr)

        patientOut.write(','+str(error))
        
    patientOut.write('\n')

patientOut.close()
print eachPatient, 'finished\n'
endTime = time.time()
workTime = endTime - startTime
print eachPatient, datetime.timedelta(seconds=workTime)

