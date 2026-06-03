# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 20:03:16 2025

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

def xInt(lineLower, lineUpper):
    offsetUpper = lineUpper.offset
    slopeUpper = lineUpper.slope

    offsetLower = lineLower.offset
    slopeLower = lineLower.slope

    voltage = (offsetUpper - offsetLower) / (slopeLower - slopeUpper)
    
    return voltage
    

#Main
voltage_s = []

B.pl.clf()
mp = B.get_file('YellowWavelengthData.txt')

voltage = mp['voltage']
current = mp['current']

B.plot_exp(x= voltage, y = current)


B.pl.ylabel('Current (A)', fontsize= 16)
B.pl.xlabel("Voltage (V)", fontsize = 16)
B.pl.title("Current vs Voltage for Yellow Light", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)

ihUpper = (0.1 <= current) & (current <= 1.5)
ihLower = (-0.004 <= current) & (current <= 0.3)

fitLower = B.linefit(voltage[ihLower], current[ihLower])
fitUpper = B.linefit(voltage[ihUpper], current[ihUpper])


voltage_s.append(xInt(fitLower, fitUpper))


#%%
B.pl.clf()
mp = B.get_file('GreenWavelengthData.txt')

voltage = mp['voltage']
current = mp['current']

B.plot_exp(x= voltage, y = current)
B.pl.ylabel('Current (A)', fontsize= 16)
B.pl.xlabel("Voltage (V)", fontsize = 16)
B.pl.title("Current vs Voltage for Green Light", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)

ihUpper = (-0.167 <= current) & (current <= 0.4)
ihLower = (-0.19 <= current) & (current <= -0.1)

fitLower = B.linefit(voltage[ihLower], current[ihLower])
fitUpper = B.linefit(voltage[ihUpper], current[ihUpper])

voltage_s.append(xInt(fitLower, fitUpper))

#%%
B.pl.clf()
mp = B.get_file('BlueWavelengthData.txt')

voltage = mp['voltage']
current = mp['current']

B.plot_exp(x= voltage, y = current)
B.pl.ylabel('Current (A)', fontsize= 16)
B.pl.xlabel("Voltage (V)", fontsize = 16)
B.pl.title("Current vs Voltage for Blue Light", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)


vhUpper = B.in_between(0.9, 1.2, voltage)
vhLower = B.in_between(1.1, 1.5, voltage)


fitLower = B.linefit(voltage[vhLower], current[vhLower])
fitUpper = B.linefit(voltage[vhUpper], current[vhUpper])

voltage_s.append(xInt(fitLower, fitUpper))

#%%
B.pl.clf()
mp = B.get_file('VioletWavelengthData.txt')

voltage = mp['voltage']
current = mp['current']

B.plot_exp(x= voltage, y = current)
B.pl.ylabel('Current (A)', fontsize= 16)
B.pl.xlabel("Voltage (V)", fontsize = 16)
B.pl.title("Current vs Voltage for Violet Light", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)

vhUpper = B.in_between(1.02, 1.4, voltage)
vhLower = B.in_between(1.2, 1.6, voltage)


fitLower = B.linefit(voltage[vhLower], current[vhLower])
fitUpper = B.linefit(voltage[vhUpper], current[vhUpper])

voltage_s.append(xInt(fitLower, fitUpper))


#%%
B.pl.clf()
mp = B.get_file('UVWavelengthData.txt')

voltage = mp['voltage']
current = mp['current']

B.plot_exp(x= voltage, y = current)
B.pl.ylabel('Current (A)', fontsize= 16)
B.pl.xlabel("Voltage (V)", fontsize = 16)
B.pl.title("Current vs Voltage for UV Light", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)

#vhUpper = B.in_between(1.2, 1.6, voltage)
#vhLower = B.in_between(1.5, 2, voltage)
irLow = B.in_between(-0.9, -0.1, current)
fitLower = B.linefit(voltage[irLow], current[irLow])

irUp = B.in_between(-0.8, 0.5, current)
fitUpper = B.linefit(voltage[irUp], current[irUp])

voltage_s.append(xInt(fitLower, fitUpper))

print(voltage_s)


#%%
B.pl.clf()
wavelength = np.array([578, 546, 436, 405, 365])
invWave = 1 / wavelength

eV = np.array(voltage_s)

B.plot_exp(x = invWave, y = eV)
B.linefit(x = invWave, y = eV)
B.pl.ylabel('Stopping Voltage (eV)', fontsize= 16)
B.pl.xlabel("1/ lambda (nm)", fontsize = 16)
B.pl.title("eV vs 1/ lambda", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)


print(voltage_s)

