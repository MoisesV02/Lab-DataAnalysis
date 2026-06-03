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
    #Puts all the diameters as an np array in meters
    diameter = fileObj['diameter'] / 100.0
    
    #Finds all the uncertainties in meters
    uncertaintyD = [(0.1 / 2) / 100.0] * len(diameter)
    uncertaintyD = np.array(uncertaintyD)
    
    #Converting voltage to kV
    voltage = fileObj['voltage']
    uncertaintyV = [(1.0 / 2.0)] * len(voltage)
    uncertaintyV = np.array(uncertaintyV)
    
    return voltage, diameter, uncertaintyD, uncertaintyV

def angle(diameter):
    R = 6.25 / 100.0
    L = 13.7 / 100.0
    
    l1 = L - R
    theta = []
    for D in diameter:
        radius = D / 2
        l2 = np.sqrt((R ** 2) - (radius ** 2))
        
        tan2theta = radius / (l1 + l2)
        theta.append((np.atan(tan2theta)) / 2.0)
        
    #Returns theta in radians
    return theta

#IDK if this is right
def uncertaintyAngle(diameter, uD):
    R = 6.25 / 100.0
    L = 13.7 / 100.0
    
    l1 = L - R
    
    uncertainty = []
    for D, u in zip(diameter, uD):
        radius = D / 2.0
        
        firstTerm = l1 + (((R**2) - (radius ** 2)) ** (1.0 / 2.0))
        secondTerm = (radius ** 2) * (1.0 / (((R**2) - (radius ** 2)) ** (1.0 / 2.0)))
        denom = (l1 + (((R**2) - (radius ** 2)) ** (1.0 / 2.0))) ** 2
        
        
        quotient = (1.0 / 4.0) * ((firstTerm - secondTerm) / denom)
        arcTan = 1.0 / (1 + (radius / (l1 + (((R ** 2) - (radius ** 2)) ** (1.0 / 2.0)))))
        
        product = quotient * arcTan
        
        uncertainty.append(np.sqrt((product * u) ** 2))
        
        
    return uncertainty

def graphUncertainties(angle, uTheta, voltage, uV):
    uncertaintyY = []
    uncertaintyX = []
    
    for theta, uY, v, uX in zip(angle, uTheta, voltage, uV):
        yPrime = np.cos(theta)
        xPrime = (1.0 / (2.0 * (v ** (3.0 / 2.0))))
        
        uncertaintyY.append(np.sqrt((yPrime * uY) ** 2))
        uncertaintyX.append(np.sqrt((xPrime * uX) ** 2))
        
    return uncertaintyY, uncertaintyX
    
    
B.pl.clf()   
voltage, diameterInner, uDiaInner, uV = diameterandUncertainty('InnerRingMeasurement.txt')
voltage, diameterOuter, uDiaOuter, uV = diameterandUncertainty('OuterRingMeasurement.txt')

print(uDiaInner)
print(uV)
print(uDiaOuter)
print(uV)


angleOuter = angle(diameterOuter)
angleInner = angle(diameterInner)
print("Outer")
print(angleOuter)
print("Inner")
print(angleInner)


angleUncertaintyInner = uncertaintyAngle(diameterInner, uDiaInner)
angleUncertaintyOuter = uncertaintyAngle(diameterOuter, uDiaOuter)



sinradInner = []
sinradOuter = []
inverseVolt = []

for v, radInner, radOuter in zip(voltage, angleInner, angleOuter):
    sinradInner.append(np.sin(radInner))
    sinradOuter.append(np.sin(radOuter))
    inverseVolt.append(1.0 / np.sqrt(v))
    

uncertaintySinInner, uncertaintyInvVolt = graphUncertainties(angleInner, angleUncertaintyInner, voltage, uV)
uncertaintySinOuter, uncertaintyInvVolt = graphUncertainties(angleOuter, angleUncertaintyOuter, voltage, uV)

print("asdfafda")
print(sinradOuter)

inverseVolt = np.array(inverseVolt)
uncertaintyInvVolt = np.array(uncertaintyInvVolt)
sinradInner = np.array(sinradInner)
sinradOuter = np.array(sinradOuter)
uncertaintySinInner = np.array(uncertaintySinInner)
uncertaintySinOuter = np.array(uncertaintySinOuter)
#%%
B.pl.clf()

#f1 = lambda x:np.cos(x)
#f2 = lambda x:np.cos(2.*x)


#Plot for inner angle
B.plot_exp(x = inverseVolt, y =sinradInner, dy = uncertaintySinInner, xerr = uncertaintyInvVolt)
#fit = B.gen_linfit([f1,f2], x = inverseVolt, y =sinradInner, yerr = uncertaintySinInner)

#%%
#Plot for outer angle
B.plot_exp(x = inverseVolt, y =sinradOuter, dy = uncertaintySinOuter, xerr = uncertaintyInvVolt)

B.pl.ylabel("sinθ")
B.pl.xlabel("Inverse of the square root of a Volt (V), " + r"$\frac{1}{\sqrt{V_a}}$")

#B.plot_line(x = inverseVolt, y =sinradOuter)

dIn = []
dOut = []
for trigInner, trigOuter, voltInv in zip(sinradInner, sinradOuter, inverseVolt):
    term1 = 1.228 / 2.0
    
    dIn.append(term1 * voltInv * (1.0 / trigInner))
    dOut.append(term1 * voltInv * (1.0 / trigOuter))
    

dIn = np.array(dIn)
dOut = np.array(dOut)
print()
print("distance with inner angle:")
print(dIn)
print("mean: ")
print( np.mean(dIn))

print()
print("Distance with outer angle:")
print(dOut)
print("mean: ")
print(np.mean(dOut))

dIn = (dIn + dOut) / 2.0

print()
print("distances summed:")
print(dIn)
print("mean: ")
print( np.mean(dIn))


    
