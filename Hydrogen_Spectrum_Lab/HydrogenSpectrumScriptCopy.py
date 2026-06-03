# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 10:52:22 2025

@author: varga
"""

import numpy as np
import LT.box as B

def convertRadians(degrees, arcMinute):
    angle = []
    
    if (isinstance(degrees, list) and isinstance(arcMinute, list)):
        for d, a in zip(degrees, arcMinute):
            theta = d + (a/60.0)
            theta = theta * (np.pi) / 180.0
            angle.append(theta)
            
    else:
        theta = degrees + (arcMinute / 60.0)
        theta = theta * (np.pi) / 180.0
        angle = theta
        
    return angle

def getData(fileName):
    fileObj = B.get_file(fileName)
    
    degrees = fileObj['DegreeFirstOrder']
    arcMin = fileObj['ArcMinuteFirstOrder']
    degreesSecondOrder = fileObj['DegreeSecondOrder']
    arcMinSecondOrder = fileObj['ArcMinuteSecondOrder']
    
    
    anglesFirstOrder = convertRadians(degrees, arcMin)
    anglesSecondOrder = convertRadians(degreesSecondOrder, arcMinSecondOrder)
    
    return anglesFirstOrder, anglesSecondOrder

def thetaOut(angle, THETA_IN):
    theta = []
    
    thetaAxis = convertRadians(198, 1)
    secondTerm = thetaAxis + THETA_IN
    
    for a in angle:
        theta.append(a - secondTerm)
        
    return np.array(theta)

def deltaTotal(GRATING_CONSTANT, THETA_IN, thetaOut):
    delta = []
    
    for a in thetaOut:
        product = GRATING_CONSTANT * (np.cos(THETA_IN) - np.cos(a))
        delta.append(product)
        
    return np.array(delta)

def wavelength(m, delta):
    wavelength = []
    for d in delta:
        #wavelength.append(d / m * (10 ** 9))
        wavelength.append(d / m)
        
    return np.array(wavelength)



#Main
#Constant in Meters
GRATING_CONSTANT = (1/1200.0) / (10 ** 3)


#Uncertainty in radians
UNCERTAINTY_VERNIER = 0.0044
UNCERTAINTY_WIDTH = 0.0025
UNCERTAINTY_ANGLE = np.sqrt((UNCERTAINTY_VERNIER ** 2) + (UNCERTAINTY_WIDTH ** 2))


THETA_IN = (convertRadians(161, 1) - convertRadians(198, 1)) / 2.0


thetaFirstOrder, thetaSecondOrder = getData('SpectralLineData.txt')

thetaOutFirst = thetaOut(thetaFirstOrder, THETA_IN)
thetaOutSecond = thetaOut(thetaSecondOrder, THETA_IN)

deltaFirst = deltaTotal(GRATING_CONSTANT, THETA_IN, thetaOutFirst)
deltaSecond = deltaTotal(GRATING_CONSTANT, THETA_IN, thetaOutSecond)

wavelengthFirst = wavelength(1.0, deltaFirst)
wavelengthSecond = wavelength(2.0, deltaSecond)

print("Wavelength first order")
print(wavelengthFirst)
print()
print("Wavelength Second order")
print(wavelengthSecond)

inverseWaveFirst = 1.0 / wavelengthFirst
inverseWaveFirst = np.array(sorted(list(inverseWaveFirst)))
print()
print("INverse Wave First")
print(inverseWaveFirst)
print( 1.0 / inverseWaveFirst)


#%%
#B.pl.clf()
nOne = 1
nTwo = np.array([nOne + 1, nOne + 2, nOne + 3, nOne + 4])
print(nTwo)
xAxis = (1 / (nOne ** 2)) - (1.0 / (nTwo ** 2))
print()
print(xAxis)
#xAxis = (1.0 / (nTwo ** 2))

xAxis = np.array(sorted(list(xAxis)))

B.plot_exp(x = xAxis, y = inverseWaveFirst)
B.linefit(x = xAxis, y = inverseWaveFirst)


#%%
inverseWaveSecond = 1.0 / wavelengthSecond
inverseWaveSecond = np.array(sorted(list(inverseWaveSecond)))



#B.pl.clf()
nOne = 2
nTwo = np.array([nOne + 1, nOne + 2, nOne + 3, nOne + 4])
print(nTwo)
xAxis = (1 / (nOne ** 2)) - (1.0 / (nTwo ** 2))
print()
print(xAxis)
#xAxis = (1.0 / (nTwo ** 2))

xAxis = np.array(sorted(list(xAxis)))

B.plot_exp(x = xAxis, y = inverseWaveSecond)
B.linefit(x = xAxis, y = inverseWaveSecond)
