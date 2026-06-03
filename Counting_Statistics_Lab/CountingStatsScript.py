# -*- coding: utf-8 -*-
#Used this for the lab report
"""
Created on Tue Sep 16 19:39:35 2025

@author: varga
"""

import numpy as np
from scipy.special import factorial
import statistics
import LT.box as B


def set_range(first_bin = 0, bin_width = 1., Nbins =10):
 """
 helper function to set the range and bin width
 input : first_bin = bin_center of the first bin,
 bin_width = the bin width,
 Nbins = total number of bins
 returns: a tuple that you can use in the range key word when defining a histogram.
 NOTE: for the histogram use the same number of bins:
 example: h = histo( r, range = set_range(-5., 1, 11), bins = 11)
this created a histogram where the first bin is centered at-5. ,
 the next at-4. etc. a total of 11 bins are created and the bin center
 of the last one is at 5. = first_bin + (Nbins-1)*bin_width
 """
 rmin = first_bin- bin_width/2.
 rmax = rmin + Nbins*bin_width
 return (rmin,rmax)

        
def poisson(array, mean):
    uniqueValueArray = set(array)
    probabilities = []
    for value in uniqueValueArray:
        numerator = (mean ** value) * (np.exp(1) ** (mean * -1))
        denominator = factorial(value)
        probabilities.append(numerator / denominator)
        
    return uniqueValueArray, probabilities

        
def gaussian(array, mean, variance):
    uniqueValueArray = set(array)
    probabilities = []
    for value in uniqueValueArray:
        exponent = ((value - mean) ** 2) / (2 * variance * -1)
        denominator = np.sqrt(2 * np.pi) * np.sqrt(variance)
        product = np.exp(exponent) / denominator
        probabilities.append(product)
        
    return uniqueValueArray, probabilities

def expectedDecays(trials, backRad):
    expectedDecays = []
    for t in trials:
        if t < 49:
            #Time duration was 40 seconds for all trials before 49.
            expectedDecays.append(backRad * 40.0)
        
        elif t >= 49 and t < 53:
            #Changed time duration from 40 seconds to 30 seconds in trials 49-52.
            expectedDecays.append(backRad * 30.0)
            
        elif t >= 53 and t < 61:
            #Changed the time duration back to 40 seconds in trials 53-61.
            expectedDecays.append(backRad * 40.0)
            
        elif t >= 61 and t < 65:
            #Changed time duration from 40 seconds  to 30 seconds in trials 61-65.
            expectedDecays.append(backRad * 30.0)
            
        elif t >= 65 and t < 71:
            #Changed time duration from 30 seconds to 20 seconds in trials 65-70.
            expectedDecays.append(backRad * 20.0)
            
        else:
            #Changed time duration from 20 seconds to 10 seconds in trial 71 onward.
            expectedDecays.append(backRad * 10.0)
    
    return expectedDecays

def uncertaintyExpectedDecays(expectedDecays):
    uncertainty = []
    for i in expectedDecays:
        uncertainty.append(np.sqrt(i))
        
    return uncertainty

def uncertaintyTotal(backRadUncertainty, uncertaintyRuler, uncertaintyExpected):
    uncertainty = []
    for u in uncertaintyExpected:
        uncertainty.append(np.sqrt((backRadUncertainty ** 2) + (u ** 2) +(uncertaintyRuler ** 2)))
        
    return uncertainty

def decaysPerSec(trials, adjustedDecays, uncertaintyTotal):
    decaysPerSec = []
    uncertainty = []
    for t, count, u in zip(trials, adjustedDecays, uncertaintyTotal):
        if t < 49:
            #Time duration was 40 seconds for all trials before 49.
            decaysPerSec.append(count / 40.0)
            uncertainty.append(u / 40.0)
        
        elif t >= 49 and t < 53:
            #Changed time duration from 40 seconds to 30 seconds in trials 49-52.
            decaysPerSec.append(count / 30.0)
            uncertainty.append(u / 30.0)
            
        elif t >= 53 and t < 61:
            #Changed time duration back to 40 seconds in trials 53-65.
            decaysPerSec.append(count / 40.0)
            uncertainty.append(u / 40.0)
            
        elif t >= 61 and t < 65:
            #Changed time duration from 40 seconds  to 30 seconds in trials 61-65.
            decaysPerSec.append(count / 30.0)
            uncertainty.append(u / 30.0)
            
        elif t >= 65 and t < 71:
            #Changed time duration from 30 seconds to 20 seconds in trials 65-70.
            decaysPerSec.append(count / 20.0)
            uncertainty.append(u / 20.0)
            
        else:
            #Changed time duration from 20 seconds to 10 seconds in trial 71 onward.
            decaysPerSec.append(count / 10.0)
            uncertainty.append(u / 10.0)
    
    return decaysPerSec, uncertainty

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
        
    
    
    

#background decay over 500 seconds
BACKGROUND_DECAY = 224.0 / 500.0
UNCERTAINTY_BACKGROUND = np.sqrt(BACKGROUND_DECAY)
UNCERTAINTY_N = 26.0
DIAMETER = 3.5 / 100.0
UNCERTAINTY_RULER = 0.05 / 100.0

mf = B.get_file('CountNumberLow.txt')
mcf = B.get_file('CountNumberHigh.txt')
mp = B.get_file('CountAtDifDistance.txt')

lowDecays = mf['Decays']
highDecays = mcf['Decays']
trials = mp['Trials']
distanceDecays = mp['Decays']
distance = mp['Distance'] / 100.0
distance = set(distance)


B.pl.clf()

voltage = np.array([700, 720, 740, 760, 780, 800])
decays = np.array([562, 658, 676, 687, 660, 736])
dy = np.array([UNCERTAINTY_N] * len(decays))

B.plot_exp(x = voltage, y = decays, dy = dy)
B.plot_line(x = voltage, y = decays)
B.pl.ylabel("Decays")
B.pl.xlabel("Voltage (V)")
B.pl.title("Detector Plateau")


#%%

#Step 1
B.pl.clf()
histoLowCount = B.histo(lowDecays, range = set_range(0, 1., 10), bins = 10, xlabel= 'Decays', ylabel='Frequency', title='Frequency of Decays')
histoLowCount.plot()

frequencyLowCount = statistics.mode(lowDecays)
meanLowCount = statistics.mean(lowDecays)
varianceLowCount = statistics.variance(lowDecays)

print("Frequency of Low Stats Data: " + str(frequencyLowCount))
print("Mean of Low Stats Data: " + str(meanLowCount))
print("Variance of Low Stats Data: " + str(varianceLowCount))

print()
histoLowCount.fit()
print()

#%%
#Step 1 Poisson graph
B.pl.clf()
uniqueLowDecays, probabilities = poisson(lowDecays, meanLowCount)

B.plot_exp(x = list(uniqueLowDecays), y = probabilities)
B.plot_line(x = list(uniqueLowDecays), y = probabilities)
B.pl.ylabel("Probability")
B.pl.xlabel("Amount of Decays")
B.pl.title("Probability of Each Decay")

#%%

#Step 2
B.pl.clf()
histoHighCount = B.histo(highDecays, range = set_range(80, 10., 10), bins = 20, xlabel= 'Decays', ylabel='Frequency', title='Frequency of Decays')
histoHighCount.plot()

frequencyHighCount = statistics.mode(highDecays)
meanHighCount = statistics.mean(highDecays)
standardDevHighCount = np.sqrt(statistics.variance(highDecays))

print("Frequency of High Stats Data: " + str(frequencyHighCount))
print("Mean of High Stats Data: " + str(meanHighCount))
print("Standard Deviation of High Stats Data: " + str(standardDevHighCount))

histoHighCount.A.set(32)
histoHighCount.mean.set(122.5)
histoHighCount.sigma.set(11.9)
histoHighCount.plot_guess()

#%%
#Step 2 Gaussian Graph
B.pl.clf()
uniqueHighDecays, probabilitiesGaussHigh = gaussian(highDecays, meanHighCount, (varianceLowCount ** 2))
B.plot_exp(x = list(uniqueHighDecays), y = probabilitiesGaussHigh)
B.plot_line(x = list(uniqueHighDecays), y = probabilitiesGaussHigh)
B.pl.ylabel("Probability")
B.pl.xlabel("Amount of Decays")
B.pl.title("Probability of Each Decay")


#%%

#Step 3
print()
print("Background Decay: " + str(BACKGROUND_DECAY))
print("Uncertainty of Background Decay: " + str(UNCERTAINTY_BACKGROUND))
print()
#%%

#Step 4
B.pl.clf()
expectedEvents = expectedDecays(trials, BACKGROUND_DECAY)
uncertaintyExpectedEvents = uncertaintyExpectedDecays(expectedEvents)


eventsAdjusted = distanceDecays - expectedEvents


uncertaintyEventsAdjusted = uncertaintyTotal(UNCERTAINTY_BACKGROUND, UNCERTAINTY_RULER, uncertaintyExpectedEvents)


eventsPerSec, uncertaintyPerSec = decaysPerSec(trials, eventsAdjusted, uncertaintyEventsAdjusted)

eventsTrial = np.split(np.array(eventsPerSec), 4)
uncertaintyTrial = np.split(np.array(uncertaintyPerSec), 4)

weightedDecays = []
weightedUncertainty = []
inverseDistance = []


for d, n, u in zip(distance, eventsTrial, uncertaintyTrial):
    inverseDistance.append(1.0 / (d ** 2))
    mean, meanUncertainty = calcMeanandUncertainty(n, u)
    weightedDecays.append(mean)
    weightedUncertainty.append(meanUncertainty)
    
inverseDistance = sorted(inverseDistance, reverse=True)
weightedDecays = sorted(weightedDecays, reverse=True)
weightedUncertainty = sorted(weightedUncertainty, reverse=True)




B.plot_exp(x = inverseDistance, y = weightedDecays, dy = weightedUncertainty, xerr = UNCERTAINTY_RULER)

line = B.linefit(x = np.array(inverseDistance), y = np.array(weightedDecays), yerr = np.array(weightedUncertainty))
B.pl.ylabel("Counts (decays/s)")
B.pl.xlabel(r"$\text{1 / }{\text{Distance}^2}$ (m$^{-2}$)")
B.pl.title("Counts over Distance Squared")

m = line.slope
print()
print(m)
m_err = np.sqrt(line.cov[0][0])

theoreticalSlopeVal = ((3.7e10) * (DIAMETER ** 2)) /  (16)
uncertaintyTheoreticalSlope = np.sqrt(theoreticalSlopeVal)

source = m / theoreticalSlopeVal
uncertaintySource = np.sqrt(((m_err / theoreticalSlopeVal) ** 2 ) + ((m * uncertaintyTheoreticalSlope / (theoreticalSlopeVal ** 2)) ** 2))
print()
print("Source Strength: " + str(source))
print("Source Uncertainty: " + str(uncertaintySource))