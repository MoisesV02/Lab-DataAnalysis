# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 11:14:37 2025

@author: varga
"""

import numpy as np
import LT.box as B

#Function that gives the weighted mean and its uncertainty.
#Returns two floats, one for the mean and one for the uncertainty.
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

#Function that gets diameter and uncertainty
def diameterandUncertainty(fileName):
    fileObj = B.get_file(fileName)
    #Puts all the diameters as an np array in cm
    diameter = fileObj['diameter']
    
    #Finds all the uncertainties in cm
    uncertaintyD = [(0.1 / 2)] * len(diameter)
    uncertaintyD = np.array(uncertaintyD)
    
    #Returns voltage
    voltage = fileObj['voltage']
    
    #Gives uncertainty of the voltage
    uncertaintyV = [(1.0 / 2.0)] * len(voltage)
    uncertaintyV = np.array(uncertaintyV)
    
    return voltage, diameter, uncertaintyD, uncertaintyV

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


def uncertaintyAngle(diameter, uD):
    R = 6.25
    L = 13.7
    
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

def dUncertainty(uSin, uInvVolt, invVolt, sintheta):
    uncertainty = []
    
    constant = 1.228 / 2.0
    for uS, uIV, iV, sT in zip(uSin, uInvVolt, invVolt, sintheta):
        firstPart = (constant * iV) / sT
        firstTerm = (uIV / iV) ** 2
        secondTerm = (uS / sT) ** 2
        
        uncertainty.append(firstPart * np.sqrt(firstTerm + secondTerm))
        
    uncertainty = np.array(uncertainty)
    
    return uncertainty
    
    
    
#Main
voltage, diameterInner, uDiaInner, uV = diameterandUncertainty('InnerRingMeasurement.txt')
voltage, diameterOuter, uDiaOuter, uV = diameterandUncertainty('OuterRingMeasurement.txt')


angleOuter = angle(diameterOuter)
angleInner = angle(diameterInner)


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



sinradInner = np.array(sinradInner)
sinradOuter = np.array(sinradOuter)
uncertaintySinInner = np.array(uncertaintySinInner)
uncertaintySinOuter = np.array(uncertaintySinOuter)
inverseVolt = np.array(inverseVolt)
uncertaintyInvVolt = np.array(uncertaintyInvVolt)

#%%
B.pl.clf()
B.plot_exp(x = inverseVolt, y = sinradInner, dy = uncertaintySinInner, xerr = uncertaintyInvVolt)
B.plot_line(x = inverseVolt, y = sinradInner)

B.plot_exp(x = inverseVolt, y = sinradOuter, dy = uncertaintySinOuter, xerr = uncertaintyInvVolt)
B.plot_line(x = inverseVolt, y = sinradOuter)

B.pl.ylabel("sinθ")
B.pl.xlabel("Inverse of the square root of a Volt (V), " + r"$\frac{1}{\sqrt{V_a}}$")
B.pl.title("sinθ vs " + r"${V_a}^{-\frac{1}{2}}$")

#%%
B.pl.clf()
sinradAvg = (sinradInner + sinradOuter) / 2.0
uncertaintySinAvg = (uncertaintySinInner + uncertaintySinOuter) / 2.0

B.plot_exp(x = inverseVolt, y = sinradAvg, dy = uncertaintySinAvg, xerr = uncertaintyInvVolt)
B.plot_line(x = inverseVolt, y = sinradAvg)



B.pl.ylabel("sinθ")
B.pl.xlabel("Inverse of the square root of a Volt (V), " + r"$\frac{1}{\sqrt{V_a}}$")
B.pl.title("Average sinθ of Inner and Outer angles  vs " + r"${V_a}^{-\frac{1}{2}}$")


#%%
dUncertainty = dUncertainty(uncertaintySinAvg, uncertaintyInvVolt, inverseVolt, sinradAvg)

d = []
for trigAvg, voltInv in zip(sinradAvg, inverseVolt):
    term1 = 1.228 / 2.0
    
    d.append(term1 * voltInv * (1.0 / (trigAvg)))
    
dMean, dU = calcMeanandUncertainty(d, dUncertainty)

print("Mean Distance between Atoms: ")
print(dMean)
print()
print("Mean Uncertainty of Distance:")
print(dU)