# -*- coding: utf-8 -*-
"""
Created on Wed Sep 24 15:28:35 2025

@author: varga
"""
import numpy as np
from scipy.special import factorial
import statistics
import LT.box as B

def expectedBDecays(trials, BACKGROUND_DECAY):
    events = []
    for t in trials:
        if t < 49:
            events.append(BACKGROUND_DECAY * 40.0)
            
        elif t >= 49 and t < 53:
            events.append(BACKGROUND_DECAY * 30.0)
            
            #Change this line if no work
        elif t >= 53 and t < 61:
            events.append(BACKGROUND_DECAY * 40.0)
            
        elif t >= 61 and t < 65:
            events.append(BACKGROUND_DECAY * 30.0)
            
        elif t >= 65 and t < 71:
            events.append(BACKGROUND_DECAY * 20.0)
            
        else:
            events.append(BACKGROUND_DECAY * 10.0)
            
    return events

def expectedUncertainty(array):
    uncertainty = []
    for i in array:
        uncertainty.append(np.sqrt(i))
        
    return uncertainty

def totalUncertainty(UNCERTAINTY_BACKGROUND, expectedUncertainty):
    uncertainty = []
    for u in expectedUncertainty:
        uncertainty.append(np.sqrt((u ** 2) + (UNCERTAINTY_BACKGROUND ** 2)))
            
    return uncertainty

def countPerSec(trials, totalDecays):
    events = []
    for t, d in zip(trials, totalDecays):
        if t < 49:
            events.append(d / 40.0)
            
        elif t >= 49 and t < 53:
            events.append(d / 30.0)
            
            #Change this line if no work
        elif t >= 53 and t < 61:
            events.append(d / 40.0)
            
        elif t >= 61 and t < 65:
            events.append(d / 30.0)
            
        elif t >= 65 and t < 71:
            events.append(d / 20.0)
            
        else:
            events.append(d / 10.0)
            
    return events

def UncertaintyPerSec(trials, totalUncertainty):
    uncertainty = []
    for t, u in zip(trials, totalUncertainty):
        if t < 49:
            uncertainty.append(u / 40.0)
            
        elif t >= 49 and t < 53:
            uncertainty.append(u / 30.0)
            
            #Change this line if no work
        elif t >= 53 and t < 61:
            uncertainty.append(u / 40.0)
            
        elif t >= 61 and t < 65:
            uncertainty.append(u / 30.0)
            
        elif t >= 65 and t < 71:
            uncertainty.append(u / 20.0)
            
        else:
            uncertainty.append(u / 10.0)
            
    return uncertainty
            

mp = B.get_file('CountAtDifDistance.txt')

trials = mp['Trials']
distanceDecays = mp['Decays']
distance = mp['Distance'] / 100.0

#Step 3
#background decay over 500 seconds
BACKGROUND_DECAY = 224.0 / 500.0
DIAMETER = 3.5 / 100.0
UNCERTAINTY_BACKGROUND = np.sqrt(BACKGROUND_DECAY)


B.pl.clf()
expectedBDecays = expectedBDecays(trials, BACKGROUND_DECAY)
uncertaintyExpectedDecays = expectedUncertainty(expectedBDecays)
#print(uncertaintyExpectedDecays)

totalDecays = distanceDecays - expectedBDecays

totalUncertainty = totalUncertainty(UNCERTAINTY_BACKGROUND, uncertaintyExpectedDecays)

decaysPerSec = countPerSec(trials, totalDecays)
uncertaintyPerSec = UncertaintyPerSec(trials, totalUncertainty)

decaysPerSecInTrials = np.split(np.array(decaysPerSec), 4)
uncertaintyPerSecInTrials = np.split(np.array(uncertaintyPerSec), 4)

averageDecaysPerSec = []
averageUncertaintyPerSec = []

for n, u in zip(decaysPerSecInTrials, uncertaintyPerSecInTrials):
    averageDecaysPerSec.append(np.mean(n))
    averageUncertaintyPerSec.append(np.mean(u))

averageDecaysPerSec = np.array(averageDecaysPerSec)
averageUncertaintyPerSec = np.array(averageUncertaintyPerSec)
inverseDistance = 1.0 /(distance ** 2)
inverseDistance = sorted(set(np.array(inverseDistance)))

#%%

#Attempt to make slope negative
inverseDistance = np.sqrt(inverseDistance)
graphDistance = []
for distance in inverseDistance:
    graphDistance.append(1 / distance)
    
B.pl.clf()
B.plot_exp(x = graphDistance, y = averageDecaysPerSec, dy = averageUncertaintyPerSec)
line = B.linefit(x = np.array(graphDistance), y = np.array(averageDecaysPerSec), yerr = np.array(averageUncertaintyPerSec))
B.pl.ylabel("Counts")
B.pl.xlabel("Distance (m)")
B.pl.title("Counts over distance")

m = line.slope
m_err = np.sqrt(line.cov[0][0])

theoreticalSlopeVal = ((3.7e10) * (DIAMETER ** 2)) /  (16)
uncertaintyTheoreticalSlope = np.sqrt(theoreticalSlopeVal)

source = m / theoreticalSlopeVal
uncertaintySource = np.sqrt(((m_err / theoreticalSlopeVal) ** 2 ) + ((m * uncertaintyTheoreticalSlope / (theoreticalSlopeVal ** 2)) ** 2))
print()
print("Source Strength: " + str(source))
print("Source Uncertainty: " + str(uncertaintySource))