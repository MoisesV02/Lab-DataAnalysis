# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 11:09:37 2025

@author: varga
"""

import numpy as np
import LT.box as B

def getDataCoarse(fileName):
    fileObj = B.get_file(fileName)
    
    degrees = fileObj['Angle']
    count = fileObj['Count']
    
    return degrees, count

def getDataFine(fileName):
    fileObj = B.get_file(fileName)
    
    degrees = fileObj['Angle']
    count = fileObj['Count']
    arcMin = fileObj['ArcMinute']
        
    degrees = degrees + (arcMin / 60.0)
    
    return degrees, count

def countUncertainty(count):
    uncertainty = []
    
    for c in count:
       uncertainty.append(np.sqrt(c)) 
    
    
    return np.array(uncertainty)

degreesCoarse, countCoarse = getDataCoarse('CoarseAnglesMeasurements.txt')
degreesFine, countFine = getDataFine('FineAngleMeasurements.txt')

all_angles = np.append(degreesCoarse, degreesFine)
all_counts = np.append(countCoarse, countFine)
uncertaintyCount = countUncertainty(all_counts) / 10.0

all_counts = all_counts / 10.0


B.pl.clf()

histo = B.histo(all_counts, bins = 120)
histo.plot()


B.pl.xlabel("Angle 2θ (°)", fontsize = 16)
B.pl.ylabel("Counts per Second", fontsize = 16)
B.pl.title("Counts per Second vs Angle", fontsize = 16)
#B.pl.tick_params(axis='both', labelsize=16)


#%%
B.pl.clf()
peakAlpha = np.array(all_angles >= 31.0) & (all_angles <= 33.0)
anglesAlpha = all_angles[peakAlpha]
countsAlpha = all_counts[peakAlpha]
uncertaintyAlpha = uncertaintyCount[peakAlpha]

B.plot_exp(x = anglesAlpha, y = countsAlpha, dy = uncertaintyAlpha)

countsAlphaMean = np.mean(countsAlpha)


#%%

