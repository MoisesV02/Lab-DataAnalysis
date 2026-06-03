# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 17:52:41 2025

@author: varga
"""

import numpy as np
import LT.box as B
import math

#Functions to calculate uncertainty of magnetic field 
def calcuncertaintyB(current, coilRadius, uncertaintyCurrent, uncertaintyradius):
    uncertainty = []
    for ampere in current:
        uncertainty.append(math.sqrt(pow((((ampere * 132) / pow(coilRadius, 2)) * uncertaintyradius), 2) + pow(((132 / coilRadius) * uncertaintyCurrent), 2)))
    
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

#Constants, didn't change in the lab. Distance measurement are in cm.
COILRADIUS = 0.296 / 2.0
UNCERTAINTY_COILRADIUS = COILRADIUS / 2.0
UNCERTAINTY_ARCPATH = 0.005 / 2.0

UNCERTAINTY_CURRENT = 0.01 / 2.0
UNCERTAINTY_VOLTAGE = 1 / 2.0
UNCERTAINTY_COIL = 0.01 / 2.0

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

#variables that have certain voltages from the experiment
sel1 = 200
sel2 = 250
sel3 = 300
sel4 = 400

#Checks to see where the voltage is 200 in the Measurement.txt file
indices = np.isin(voltage, sel1)

B.pl.clf()

B.plot_exp(x = current, y = ratioEM, dy = uncertaintyRatio, xerr = UNCERTAINTY_CURRENT)
line = B.linefit(x = current, y = ratioEM, yerr = uncertaintyRatio)

#B.linefit(x = current[indices], y = ratioEM[indices], yerr = uncertaintyRatio[indices])

#%%
indices = np.isin(voltage, sel2)
#B.linefit(x = current[indices], y = ratioEM[indices], yerr = uncertaintyRatio[indices])

#%%
indices = np.isin(voltage, sel3)
#B.linefit(x = current[indices], y = ratioEM[indices], yerr = uncertaintyRatio[indices])

#%%
indices = np.isin(voltage, sel4)
#B.linefit(x = current[indices], y = ratioEM[indices], yerr = uncertaintyRatio[indices])

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

print("First row:")
print("Voltage:", voltage[0])
print("Current:", current[0])
print("Diameter (m):", diameter[0])
print("Radius (m):", radius[0])
print("B (T):", magneticFieldB[0])
print("Uncertainty B:", uncertaintyB[0])
print("e/m:", ratioEM[0])
print("Uncertainty e/m:", uncertaintyRatio[0])
