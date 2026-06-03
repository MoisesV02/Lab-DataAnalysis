# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 17:52:41 2025

@author: varga
"""

import numpy as np
import LT.box as B
import math

#Functions to calculate uncertainty of magnetic field 
def calcuncertaintyB(current, radius, uncertaintyCurrent, uncertaintyradius):
    uncertainty = []
    for ampere in current:
        uncertainty.append(math.sqrt(pow((((ampere * 132) / pow(radius, 2)) * uncertaintyradius), 2) + pow(((132 / radius) * uncertaintyCurrent), 2)))
    
    return uncertainty

#Function to calculate the ratio between charge & mass
def calcuncertaintyRatio(voltage, magneticField, radiusEPath, uncertaintyVoltage, uncertaintyField, uncertaintyPath):
    uncertainty = []
    quadV = 0
    firstBR = 0
    secondBR = 0
    thirdBR = 0
    for v, b, r, uB in zip(voltage, magneticField, radiusEPath, uncertaintyField):
        quadV = 4 * v
        firstBR = (b ** 2) * (r ** 2)
        secondBR = (b ** 3) * (r ** 2)
        thirdBR = (b ** 2) * (r ** 3)
        uncertainty.append(math.sqrt(pow(((2.0 / firstBR) * uncertaintyVoltage), 2) + pow(((quadV / secondBR) * uB), 2) + pow(((quadV / thirdBR) * uncertaintyPath), 2)))
    
    return uncertainty

#Function to calculate weighted average and uncertainty. Returns both as a tuple.
def calcMeanandUncertainty(ratioEM, uncertaintyRatio):
    weightedSum = 0.0
    weightedTotal = 0.0
    for rEM, uRatio in zip(ratioEM, uncertaintyRatio):
        weight = 1.0 / (uRatio ** 2)
        weightedSum += rEM * weight
        weightedTotal += weight
        
    mean = weightedSum / weightedTotal
    uncertainty = math.sqrt(1.0 / weightedTotal)
       
    
    return mean, uncertainty
        

        
#Main
mf = B.get_file('Measurements.txt')

#Constants, didn't change in the lab. Distance measurement are in m.
COILRADIUS = 0.296 / 2.0
UNCERTAINTY_COILRADIUS = COILRADIUS / 2.0
UNCERTAINTY_ARCPATH = 0.005 / 2.0

UNCERTAINTY_CURRENT = 0.01 / 2.0
UNCERTAINTY_VOLTAGE = 1 / 2.0

MAGNETIC_PERMEABILITY =  (4 * np.pi) * pow(10, -7)

#Data from the Measurements.txt file are turned into numpy arrays
voltage = mf['Voltage']
current = mf ['Current']
diameter = mf ['Diameter'] / 100.0

#Calculate the radius for each data point based on the diameter of each measurement.
radius = diameter / 2.0

#Calculate the magnetic field based on the given equation, then calls the uncertaintyB function to get the uncertainty into an array
magneticFieldB = ((132 * current * MAGNETIC_PERMEABILITY) / COILRADIUS) * (pow( 4.0 / 5.0, 3/2))
uncertaintyB = np.array(calcuncertaintyB(current, COILRADIUS, UNCERTAINTY_CURRENT, UNCERTAINTY_COILRADIUS))
uncertaintyB = MAGNETIC_PERMEABILITY * pow((4.0/5.0), 3/2) * uncertaintyB

#Calculates the charge to mass ratio of the electrons and then calls the uncertaintyRatio to get its respective uncertainty.
ratioEM = (2 * voltage) / (pow(magneticFieldB, 2) * pow(radius, 2))
uncertaintyRatio = np.array(calcuncertaintyRatio(voltage, magneticFieldB, radius, UNCERTAINTY_VOLTAGE, uncertaintyB, UNCERTAINTY_ARCPATH))


B.pl.clf()

B.plot_exp(x = current, y = ratioEM, dy = uncertaintyRatio, xerr = UNCERTAINTY_CURRENT)
line = B.linefit(x = current, y = ratioEM, yerr = uncertaintyRatio)

B.pl.xlabel("Current (A)")
B.pl.ylabel("Ratio (C/kg)")
B.pl.title("Charge to Mass vs Current")

#%%
weightedAvg, avgUncertainty = calcMeanandUncertainty(ratioEM, uncertaintyRatio)

print()
print("Magnetic field (B) array:")
print(magneticFieldB)

print()
print("Magnetic field (B) uncertainty array:")
print(uncertaintyB)

print()
print("EM Ratio:")
print(ratioEM)

print()
print("EM Ratio Uncertainty:")
print(uncertaintyRatio)

print()
print("Weighted Mean of Ratio: " + f"{weightedAvg:E}")
print("Uncertainty of Weighted Mean: " + f"{avgUncertainty:E}")

print()
print("Chi Square of Plot: " + str(line.chi_red))