# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 11:16:29 2025

@author: varga
"""

import numpy as np
import LT.box as B

def calcMeanandUncertainty(array, uncertainty):
    weightedSum = 0.0
    weightedTotal = 0.0
    for n, u in zip(array, uncertainty):
        weight = 1.0 / (u ** 2)
        weightedSum += n * weight
        weightedTotal += weight
        
    mean = weightedSum / weightedTotal
    uncertainty = np.sqrt(1.0 / weightedTotal)
       
    
    return mean, uncertainty


 # get the live time from a spectrum
def get_time(sp):
    tt = float(sp.title.split('live time:')[1].split('s,')[0].strip())
    return tt

#Main
B.pl.clf()
cal = B.get_spectrum("0_deg_empty.Spe")
speGold0 = B.get_spectrum("0_deg.Spe")

#Orange speGold0 and blue cal
cal.plot()
speGold0.plot()

#%%
# Step 1
FILES = ["5_deg.Spe", "n5_deg.Spe", "10_deg.Spe", "n10_deg.Spe", 
         "15_deg.Spe", "n15_deg.Spe", "20_deg.Spe", "n20_deg.Spe", 
         "n30_deg.Spe", "n40_deg.Spe"] #Slow rate



