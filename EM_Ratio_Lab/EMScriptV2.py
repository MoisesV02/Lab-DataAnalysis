# -*- coding: utf-8 -*-
"""
Created on Sat Sep 13 17:23:33 2025

@author: varga
"""

import numpy as np
import LT.box as B

def uncertaintyB(current, coilRadius, uncertaintyRadius, uncertaintyCurrent, constant, turns):
    uncertainty = []
    factor = constant * ((4.0 / 5.0) ** (3.0 / 2.0))
    for ampere in current:
        term1 = ((uncertaintyRadius * ampere * turns) / (coilRadius ** 2)) ** 2
        term2 = ((uncertaintyCurrent * turns) / coilRadius) ** 2
        uncertainty.append(factor * np.sqrt(term1 + term2))
        
    return uncertainty

def uncertaintyRatio(field, arcPath, uncertaintyV, voltage, uncertaintyB, uncertaintyArcPath):
    uncertainty = []
    for b, p, v, uB in zip(field, arcPath, voltage, uncertaintyB):
        term1 = ((2 * uncertaintyV) / ((b ** 2) * (p ** 2))) ** 2
        term2 = ((4 * v * uB) / ((b ** 3) * (p ** 2))) ** 2
        term3 = ((4 * v * uncertaintyArcPath) / ((b ** 2) * (p ** 3))) ** 2
        
        uncertainty.append(np.sqrt(term1 + term2 + term3))
        
    return uncertainty

def calcMeanandUncertainty(ratioEM, uncertaintyRatio):
    weightedSum = 0.0
    weightedTotal = 0.0
    for rEM, uRatio in zip(ratioEM, uncertaintyRatio):
        weight = 1.0 / (uRatio ** 2)
        weightedSum += rEM * weight
        weightedTotal += weight
        
    mean = weightedSum / weightedTotal
    uncertainty = np.sqrt(1.0 / weightedTotal)
       
    
    return mean, uncertainty


mf = B.get_file('Measurements.txt')

voltage = mf['Voltage']
current = mf['Current']

#In the file it is measured in cm
diameter = mf['Diameter'] / 100.0

radiusPath = diameter / 2.0

COILS = 132

#Originally in cm, so changed to m
COILRADIUS = 0.296 / 2.0
UNCERTAINTY_COILRADIUS = COILRADIUS / 2.0
UNCERTAINTY_ARCPATH = 0.005 / 2.0

UNCERTAINTY_CURRENT = 0.01 / 2.0
UNCERTAINTY_VOLTAGE = 1.0 / 2.0

MAGNETIC_PERMEABILITY =  (4 * np.pi) * pow(10, -7)


magneticFieldB = ((4.0 / 5.0) ** (3.0 / 2.0)) * ((COILS * current * MAGNETIC_PERMEABILITY) / COILRADIUS)
uncertaintyB = np.array(uncertaintyB(current, COILRADIUS,UNCERTAINTY_COILRADIUS, UNCERTAINTY_CURRENT, MAGNETIC_PERMEABILITY, COILS))


print(magneticFieldB)
print(uncertaintyB)

ratioEM = (2 * voltage) / ((magneticFieldB ** 2) * (radiusPath ** 2))
print(ratioEM)

uncertaintyRatio = np.array(uncertaintyRatio(magneticFieldB, radiusPath, UNCERTAINTY_VOLTAGE, voltage, uncertaintyB, UNCERTAINTY_ARCPATH))
print()
print(uncertaintyRatio)

B.pl.clf()
B.plot_exp(current, ratioEM, dy = uncertaintyRatio, xerr = UNCERTAINTY_CURRENT)
B.linefit(current, ratioEM, yerr = uncertaintyRatio)

weightedAvg, avgUncertainty = calcMeanandUncertainty(ratioEM, uncertaintyRatio)
print("Weighted Mean of Ratio: " + f"{weightedAvg:E}")
print("Uncertainty of Weighted Mean: " + f"{avgUncertainty:E}")

