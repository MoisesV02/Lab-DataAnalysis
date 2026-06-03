# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 20:03:16 2025

@author: varga
"""

import numpy as np
import LT.box as B

def errorV(oU, oL, mU, mL, uOU, uOL, uMU, uML):
    term1 = (uOU / (mL - mU)) ** 2
    term2 = (uOL / (mL - mU)) ** 2
    term3 = (uMU * (oU - oL) / ((mL - mU) ** 2)) ** 2
    term4 = (uML * (oU - oL) / ((mL - mU) ** 2)) ** 2
    
    uncertainty = np.sqrt(term1 + term2 + term3 + term4)
    
    return uncertainty
    

def xInt(lineLower, lineUpper):
    offsetUpper = lineUpper.offset
    uOffsetUpper = lineUpper.sigma_o
    slopeUpper = lineUpper.slope
    uSlopeUpper = lineUpper.sigma_s

    offsetLower = lineLower.offset
    uOffsetLower = lineLower.sigma_o
    slopeLower = lineLower.slope
    uSlopeLower = lineLower.sigma_s

    voltage = (offsetUpper - offsetLower) / (slopeLower - slopeUpper)
    error = errorV(offsetUpper, offsetLower, slopeUpper, slopeLower, uOffsetUpper, uOffsetLower, uSlopeUpper, uSlopeLower)
    
    return voltage, error
    

#Main
voltage_s = []
uVoltage_s = []

UNCERTAINTY_VOLTAGE = UNCERTAINTY_CURRENT = 5e-6

B.pl.clf()
mp = B.get_file('YellowWavelengthData.txt')

voltage = mp['voltage']
current = mp['current']

B.plot_exp(x= voltage, y = current, dy= [UNCERTAINTY_CURRENT] * len(current), xerr = [UNCERTAINTY_VOLTAGE] * len(voltage))


B.pl.ylabel('Current (mA)', fontsize= 16)
B.pl.xlabel("Voltage (V)", fontsize = 16)
B.pl.title("Current vs Voltage for Yellow Light", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)

ihUpper = (0.1 <= current) & (current <= 1.5)
ihLower = (-0.004 <= current) & (current <= 0.3)

fitLower = B.linefit(voltage[ihLower], current[ihLower])
fitUpper = B.linefit(voltage[ihUpper], current[ihUpper])

stoppingV, uV = xInt(fitLower, fitUpper)
voltage_s.append(stoppingV)
uVoltage_s.append(uV)


#%%
B.pl.clf()
mp = B.get_file('GreenWavelengthData.txt')

voltage = mp['voltage']
current = mp['current']

B.plot_exp(x= voltage, y = current, dy= [UNCERTAINTY_CURRENT] * len(current), xerr = [UNCERTAINTY_VOLTAGE] * len(voltage))
B.pl.ylabel('Current (mA)', fontsize= 16)
B.pl.xlabel("Voltage (V)", fontsize = 16)
B.pl.title("Current vs Voltage for Green Light", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)

vhUpper = (0.6 <= voltage) & (voltage <= 0.85)
vhLower = (0.7 <= voltage) & (voltage <= 1.4)


fitLower = B.linefit(voltage[vhLower], current[vhLower])
fitUpper = B.linefit(voltage[vhUpper], current[vhUpper])

stoppingV, uV = xInt(fitLower, fitUpper)

voltage_s.append(stoppingV)
uVoltage_s.append(uV)

#%%
B.pl.clf()
mp = B.get_file('BlueWavelengthData.txt')

voltage = mp['voltage']
current = mp['current']

B.plot_exp(x= voltage, y = current, dy= [UNCERTAINTY_CURRENT] * len(current), xerr = [UNCERTAINTY_VOLTAGE] * len(voltage))
B.pl.ylabel('Current (mA)', fontsize= 16)
B.pl.xlabel("Voltage (V)", fontsize = 16)
B.pl.title("Current vs Voltage for Blue Light", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)


vhUpper = B.in_between(0.9, 1.2, voltage)
vhLower = B.in_between(1.06, 1.41, voltage)


fitLower = B.linefit(voltage[vhLower], current[vhLower])
fitUpper = B.linefit(voltage[vhUpper], current[vhUpper])

stoppingV, uV = xInt(fitLower, fitUpper)
voltage_s.append(stoppingV)
uVoltage_s.append(uV)

#%%
B.pl.clf()
mp = B.get_file('VioletWavelengthData.txt')

voltage = mp['voltage']
current = mp['current']

B.plot_exp(x= voltage, y = current, dy= [UNCERTAINTY_CURRENT] * len(current), xerr = [UNCERTAINTY_VOLTAGE] * len(voltage))
B.pl.ylabel('Current (mA)', fontsize= 16)
B.pl.xlabel("Voltage (V)", fontsize = 16)
B.pl.title("Current vs Voltage for Violet Light", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)

vhUpper = B.in_between(1.02, 1.3, voltage)
vhLower = B.in_between(1.2, 1.45, voltage)


fitLower = B.linefit(voltage[vhLower], current[vhLower])
fitUpper = B.linefit(voltage[vhUpper], current[vhUpper])

stoppingV, uV = xInt(fitLower, fitUpper)
voltage_s.append(stoppingV)
uVoltage_s.append(uV)


#%%
B.pl.clf()
mp = B.get_file('UVWavelengthData.txt')

voltage = mp['voltage']
current = mp['current']

B.plot_exp(x= voltage, y = current, dy= [UNCERTAINTY_CURRENT] * len(current), xerr = [UNCERTAINTY_VOLTAGE] * len(voltage))
B.pl.ylabel('Current (mA)', fontsize= 16)
B.pl.xlabel("Voltage (V)", fontsize = 16)
B.pl.title("Current vs Voltage for UV Light", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)

vhUpper = B.in_between(1.2, 1.6, voltage)
vhLower = B.in_between(1.5, 2, voltage)


fitLower = B.linefit(voltage[vhLower], current[vhLower])
fitUpper = B.linefit(voltage[vhUpper], current[vhUpper])

stoppingV, uV = xInt(fitLower, fitUpper)
voltage_s.append(stoppingV)
uVoltage_s.append(uV)


#%%
B.pl.clf()
print()
wavelength = np.array([578, 546, 436, 405, 365])
invWave = 1 / wavelength

eV = np.array(voltage_s)
ueV = np.array(uVoltage_s)


B.plot_exp(x = invWave, y = eV, dy=ueV)
line = B.linefit(x = invWave, y = eV, yerr = ueV)
B.pl.ylabel('Stopping Voltage (eV)', fontsize= 16)
B.pl.xlabel(r"1/$\lambda$ (nm)", fontsize=16)
B.pl.title(r"eV vs 1/$\lambda$", fontsize=16)
B.pl.tick_params(axis='both', labelsize=16)


print()
print("hc from slope: " + str(line.slope))
print("Slope uncertainty: " + str(line.sigma_s))
print("Work function: " + str(line.offset))
print("Work Uncertainty: " + str(line.sigma_o))
