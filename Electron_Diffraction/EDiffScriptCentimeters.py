# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 11:14:37 2025

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

def diameterandUncertainty(fileName):
    fileObj = B.get_file(fileName)
    #Puts all the diameters as an np array in cm
    diameter = fileObj['diameter']
    
    #Finds all the uncertainties in cm
    uncertaintyD = [0.1 / 2] * len(diameter)
    uncertaintyD = np.array(uncertaintyD)
    
    return diameter, uncertaintyD

def angle(diameter):
    R = 6.25
    L = 13.7
    
    l1 = L - R
    theta = []
    for D in diameter:
        radius = D / 2
        l2 = np.sqrt((R ** 2) - (radius ** 2))
        
        tan2theta = radius / (l1 + l2)
        theta.append((np.atan(tan2theta)) / 2.0)
        
    #Returns theta in radians
    return theta
    
    
    
diameterInner, uDiaInner = diameterandUncertainty('InnerRingMeasurement.txt')
diameterOuter, uDiaOuter = diameterandUncertainty('OuterRingMeasurement.txt')
print(diameterInner)
print(uDiaInner)
print()

print(diameterOuter)
print(uDiaOuter)
print()

angleOuter = angle(diameterOuter)
print(len(angleOuter))
print(angleOuter)



#sinrad = []
#inverseVolt = []




#B.plot_exp(x = inverseVolt, y =sinrad)



