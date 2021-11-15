import pandas as pd
import re, sys
import os
import numpy as np
import random
import scipy.stats as stats
from scipy.stats import norm
import pybedtools
import bx
#from bx.intervals.io import GenomicIntervalReader
from bx.bbi.bigwig_file import BigWigFile

import pyBigWig

if len(sys.argv) != 3:
    print "Usage: python bigwig2bin.py [bigwig folder] [sampleList]"
    sys.exit(1)


tempPath=os.environ['PWD']
bigwigFolder = str(sys.argv[1])
sampleInfo = pd.read_csv(str(sys.argv[2]), index_col=0)

peakLoc = pd.read_csv(tempPath+'/source/peakLocation.csv', index_col=0).iloc[:,0:3]

outBin1 = open(tempPath+'/.tmp/bin1.csv', 'w')
outBin1.write('sampleID,'+','.join(peakLoc.index.tolist())+'\n')
outBin2 = open(tempPath+'/.tmp/bin2.csv', 'w')
outBin2.write('sampleID,'+','.join(peakLoc.index.tolist())+'\n')
outBin3 = open(tempPath+'/.tmp/bin3.csv', 'w')
outBin3.write('sampleID,'+','.join(peakLoc.index.tolist())+'\n')
outBin4 = open(tempPath+'/.tmp/bin4.csv', 'w')
outBin4.write('sampleID,'+','.join(peakLoc.index.tolist())+'\n')
outBin5 = open(tempPath+'/.tmp/bin5.csv', 'w')
outBin5.write('sampleID,'+','.join(peakLoc.index.tolist())+'\n')

for sample in sampleInfo.index:
    fileName = sample+'.bigwig'
#    inputBg = open(bigwigFolder+'/'+fileName, 'r')
#    bw = BigWigFile(inputBg)
    bw=pyBigWig.open(bigwigFolder+'/'+fileName)
    
    outBin1.write(sample)
    outBin2.write(sample)
    outBin3.write(sample)
    outBin4.write(sample)
    outBin5.write(sample)
    
    for eachPk in peakLoc.index:
        #     pkBed = pybedtools.BedTool.from_dataframe(brcaPeakLoc.loc[[eachPk]])
        pkChr = peakLoc.loc[eachPk, 'seqnames']
        pkStart = peakLoc.loc[eachPk, 'start']
        pkEnd = peakLoc.loc[eachPk, 'end']
        #     print len(bw.get(pkChr, pkStart, pkEnd))
#        print pkChr, pkStart, pkEnd

#        bin1 = [x[2] for x in bw.get(pkChr, pkStart, pkStart+100)]
#        bin2 = [x[2] for x in bw.get(pkChr, pkStart+100, pkStart+200)]
#        bin3 = [x[2] for x in bw.get(pkChr, pkStart+200, pkStart+300)]
#        bin4 = [x[2] for x in bw.get(pkChr, pkStart+300, pkStart+400)]
#        bin5 = [x[2] for x in bw.get(pkChr, pkStart+400, pkEnd)]

        bin1 = np.nan_to_num(bw.values(pkChr, pkStart, pkStart+100, numpy=True), copy=True)
        bin2 = np.nan_to_num(bw.values(pkChr, pkStart+100, pkStart+200, numpy=True), copy=True)
        bin3 = np.nan_to_num(bw.values(pkChr, pkStart+200, pkStart+300, numpy=True), copy=True)
        bin4 = np.nan_to_num(bw.values(pkChr, pkStart+300, pkStart+400, numpy=True), copy=True)
        bin5 = np.nan_to_num(bw.values(pkChr, pkStart+400, pkEnd, numpy=True), copy=True)
        

        outBin1.write(','+str(np.mean(bin1)))
        outBin2.write(','+str(np.mean(bin2)))
        outBin3.write(','+str(np.mean(bin3)))
        outBin4.write(','+str(np.mean(bin4)))
        outBin5.write(','+str(np.mean(bin5)))

    outBin1.write('\n')
    outBin2.write('\n')
    outBin3.write('\n')
    outBin4.write('\n')
    outBin5.write('\n')

outBin1.close()
outBin2.close()
outBin3.close()
outBin4.close()
outBin5.close()
