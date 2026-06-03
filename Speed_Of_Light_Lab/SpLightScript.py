# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 19:31:53 2025

@author: varga
"""

import numpy as np
import LT.box as B

def getFrequencyAndShift(fileName):
    fileObj = B.get_file(fileName)
    frequency = fileObj['frequency']
    
    distance = fileObj['distance'] *(1.0/100)  / (10 ** 3)
    
    
    return np.array(frequency), np.array(distance)
    

def omega(frequency):
    omega = []
    for f in frequency:
        omega.append(2.0 * (np.pi) * f)
        
    return np.array(omega)

def calCAndUncertainty(omega, shift, uOmega, uShift):
    #Will take the combine data set which is half negative and have positive
    #So will take absolute value, since each measurement is independent from the next and later average them
    c = []
    uncertainty = []
    omega = np.abs(omega)
    
    F = 0.252
    D1 = 24.39 / 3.281
    D2 = 23.88 / 3.281
    
    DTotal = D1 + D2

    
    for w, x in zip(omega, shift):
        numerator = 4 * DTotal * F * w
    
        c.append(numerator / x)
        
    #Calculate the uncertainty for each c
    uncertaintyD = (0.01 / 2) / 3.281
    uncertaintyDTotal = np.sqrt((2 * uncertaintyD) ** 2)
    
    for w, x, uW, uX in zip (omega, shift, uOmega, uShift):
        firstProduct = ((4 * F * DTotal * uW) / x) ** 2
        secondProduct = ((4 * w * F * DTotal *uX) / (x ** 2)) ** 2
        thirdProduct = ((4 * w * F * uncertaintyDTotal) / x) ** 2
        
        uncertainty.append(np.sqrt(firstProduct + secondProduct + thirdProduct))
        
    return c, uncertainty
        
    
    
def slopeCAndUncertainty(slope, uSlope):
    F = 0.252
    D1 = 24.39 / 3.281
    D2 = 23.88 / 3.281
    
    DTotal = D1 + D2
    
    c = (4 * F * DTotal) / slope
    
    #Uncertainty
    uncertaintyD = (0.01 / 2) / 3.281
    uncertaintyDTotal = np.sqrt((2 * uncertaintyD) ** 2)
    
    firstProduct = ((4 * F * DTotal * uSlope) / (slope ** 2)) ** 2
    secondProduct = ((4 * F * uncertaintyDTotal) / slope) ** 2
    
    uncertainty = np.sqrt(firstProduct + secondProduct)
    
    return c, uncertainty
    
    

#Main
B.pl.clf()

frequencyCW, shiftCW = getFrequencyAndShift('CWDistanceMeasurements.txt')


UNCERTAINTY_SHIFT = [10 / 4.0] * len(shiftCW) 
UNCERTAINTY_SHIFT = np.array(UNCERTAINTY_SHIFT) * ((1.0/100)  / (10 ** 3))

UNCERTAINTY_FREQUENCY = [1.0 / 2.0] * len(frequencyCW)
UNCERTAINTY_FREQUENCY = np.array(UNCERTAINTY_FREQUENCY)

UNCERTAINTY_OMEGA = 2 * (np.pi) * UNCERTAINTY_FREQUENCY

#%%
#Clockwise rotation graph, assumed clockwise was positive
omegaCW = omega(frequencyCW)


B.pl.clf()

B.plot_exp(x = omegaCW, y = shiftCW, dy = UNCERTAINTY_SHIFT, xerr = UNCERTAINTY_OMEGA)
lineCW = B.linefit(x = omegaCW, y = shiftCW, yerr = UNCERTAINTY_SHIFT)

B.pl.xlabel('Angular Velocity (rad/s)', fontsize = 14)
B.pl.ylabel('Position ∆x (m)', fontsize = 14)
B.pl.title("Angular Velocity vs Shift of CW Data", fontsize = 14)
B.pl.tick_params(axis='both', labelsize=14)

print()
print("CW slope")
slopeCW = lineCW.slope
uSlopeCW = lineCW.sigma_s
print(slopeCW)

print()
print("CW slope uncertainty")
print(uSlopeCW)

cSlopeCW, uCSlopeCW = slopeCAndUncertainty(slopeCW, uSlopeCW)

print()
print("speed of light CW")
print(f"{cSlopeCW:.3e}")

print()
print("Uncertainty of the speed of light CW")
print(f"{uCSlopeCW:.3e}")
print()




#%%
B.pl.clf()

frequencyCCW, shiftCCW = getFrequencyAndShift('CCWDistanceMeasurements.txt')


#Graph for Counter clockwise, negative because assumed that clockwise was positive
omegaCCW = omega(frequencyCCW) * (-1)


B.plot_exp(x = omegaCCW, y = shiftCCW, dy = UNCERTAINTY_SHIFT, xerr = UNCERTAINTY_OMEGA)
lineCCW = B.linefit(x = omegaCCW, y = shiftCCW, yerr = UNCERTAINTY_SHIFT)

B.pl.xlabel('Angular Velocity (rad/s)', fontsize = 14)
B.pl.ylabel('Position ∆x (m)', fontsize = 14)
B.pl.title("Angular Velocity vs Shift of CCW Data", fontsize = 14)
B.pl.tick_params(axis='both', labelsize=14)

print()
slopeCCW = lineCCW.slope
uSlopeCCW = lineCCW.sigma_s

print("CCW slope")
print(slopeCCW)

print()
print("CCW slope uncertainty")
print(uSlopeCCW)
print()

cSlopeCCW, uCSlopeCCW = slopeCAndUncertainty(slopeCCW, uSlopeCCW)
print("Speed of light CCW")
print(f"{cSlopeCCW:.3e}")

print()
print("Uncertainty of speed of light CCW")
print(f"{uCSlopeCCW:.3e}")
print()


#%%
#CW + CCW Graph
B.pl.clf()

omegaCombined = list(omegaCCW) + list(omegaCW)
shiftCombined = list(shiftCCW) + list(shiftCW)
uOmegaCombined = list(UNCERTAINTY_OMEGA) * 2
uShiftCombined = list(UNCERTAINTY_SHIFT) * 2

omegaCombined = np.array(omegaCombined)
shiftCombined = np.array(shiftCombined)
uOmegaCombined = np.array(uOmegaCombined)
uShiftCombined = np.array(uShiftCombined)

B.plot_exp(x = omegaCombined, y = shiftCombined, dy = uShiftCombined, xerr = uOmegaCombined)
lineCombined = B.linefit(x = omegaCombined, y = shiftCombined, yerr = uShiftCombined)

B.pl.xlabel('Angular Velocity (rad/s)', fontsize = 14)
B.pl.ylabel('Position ∆x (m)', fontsize = 14)
B.pl.title("Angular Velocity vs Shift of CW & CCW Data", fontsize = 14)
B.pl.tick_params(axis='both', labelsize=14)
#%%
slopeCombined = lineCombined.slope
uSlopeCombined = lineCombined.sigma_s

print()
print("slope Combined")
print(slopeCombined)
print()

print("Uncertainty slope combined")
print(uSlopeCombined)
print()

cSlopeCombined, uCSlopeCombined = slopeCAndUncertainty(slopeCombined, uSlopeCombined)
print()
print("Speed of light of both graphs combined")
print(f"{cSlopeCombined:.3e}")
print()

print("Uncertainty of speed of light of both graphs combined")
print(f"{uCSlopeCombined:.3e}")

#%%
#Calculating the speed of light and uncertainty for each data point and later averaging it
cPerData, uCPerData = calCAndUncertainty(omegaCombined, shiftCombined, uOmegaCombined, uShiftCombined)

cPerData = np.array(cPerData)
uCPerData = np.array(uCPerData)

print()
print("Mean from each measurement")
print(f"{np.mean(cPerData):.3e}")

print()
print("Average Uncertainty")
print(f"{np.mean(uCPerData):.3e}")