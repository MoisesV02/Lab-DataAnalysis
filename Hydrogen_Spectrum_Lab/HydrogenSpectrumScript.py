# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 10:52:22 2025

@author: varga
"""
#The one used in the lab report was V2
import numpy as np
import LT.box as B

def convertRadians(degrees, arcMinute):
    angle = []
    
    for d, a in zip(degrees, arcMinute):
        theta = d + (a/60.0)
        theta = theta * (np.pi) / 180.0
        angle.append(theta)
        
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

    
    
def calcThetaOutAndUncertainty(angle, THETA_IN, UNCERTAINTY_ANGLE, UNCERTAINTY_IN):
    thetaOut = []
    thetaUncertainty = []
    
    thetaAxis = 198 + (1.0/60.0)
    thetaAxis = thetaAxis * (np.pi) / 180.0
    
    for a in angle:
        thetaOut.append(a - (thetaAxis + THETA_IN))
        
        product = np.sqrt((2 * (UNCERTAINTY_ANGLE ** 2)) + (UNCERTAINTY_IN ** 2))
        thetaUncertainty.append(product)
        
    return thetaOut, thetaUncertainty
    
def deltaTotalAndUncertainty(THETA_IN, GRATING_CONSTANT, thetaOut, UNCERTAINTY_IN, uOut):
    delta = (GRATING_CONSTANT) * (np.cos(THETA_IN) - np.cos(thetaOut))
    
    uncertainty = []
    for a, u in zip(thetaOut, uOut):
        product = GRATING_CONSTANT * np.sqrt(((UNCERTAINTY_IN * np.sin(THETA_IN)) ** 2) + ((u * np.sin(a)) ** 2))
        uncertainty.append(product)
        
    return np.array(delta), np.array(uncertainty)

def wavelengthAndUncertainty(m, delta, uDelta):
    wavelength = []
    uncertainty = []
    
    for d, uD in zip(delta, uDelta):
        wavelength.append((d / m))
        uncertainty.append((uD / m))
        
    return np.array(wavelength), np.array(uncertainty)

def energyAndUncertainty(wavelength, uWave):
    energy = []
    uncertainty = []
    
    h = 6.626 * (10 ** -34)
    c = 3 * (10 ** 8)
    for w, uW in zip(wavelength, uWave):
        energy.append((h * c) / w)
        uncertainty.append((h * c * uW) / (w ** 2))
        
    return np.array(energy), np.array(uncertainty)

def inverseLambdaUncertainty(inverseWave, uWave):
    uncertainty = []
    for i, u in zip(inverseWave, uWave):
        uncertainty.append(u / ((1.0 / i) ** 2))
        
    return np.array(uncertainty)

#Main
#Constant in Meters
GRATING_CONSTANT = (1/1200.0) / (10 ** 3)


#Uncertainty in radians
UNCERTAINTY_VERNIER = 0.0044
UNCERTAINTY_WIDTH = 0.0025
UNCERTAINTY_ANGLE = np.sqrt((UNCERTAINTY_VERNIER ** 2) + (UNCERTAINTY_WIDTH ** 2))


THETA_IN = ((161 + (1/60.0)) - (198 + (1/60.0))) / 2.0
THETA_IN = THETA_IN * np.pi / 180.0
UNCERTAINTY_IN = UNCERTAINTY_ANGLE / (np.sqrt(2))



thetaFirstOrder, thetaSecondOrder = getData('SpectralLineData.txt')


thetaOutFirstOrder, thetaOutFirstOrderUncertainty = calcThetaOutAndUncertainty(thetaFirstOrder, THETA_IN, UNCERTAINTY_ANGLE, UNCERTAINTY_IN)
print(thetaOutFirstOrder)
print(THETA_IN)

delta, deltaUncertainty = deltaTotalAndUncertainty(THETA_IN, GRATING_CONSTANT, thetaOutFirstOrder, UNCERTAINTY_IN, thetaOutFirstOrderUncertainty)
print("Delta total")
print(delta)

print()
print("Delta Uncertainty")
print(deltaUncertainty * (10 ** 9))

wavelengthFirstOrder, uWaveFirstOrder = wavelengthAndUncertainty(1, delta, deltaUncertainty)
print()
print("Wavelength for First Order")
print(wavelengthFirstOrder)

print()
print("Uncertainty Wavelength")
print(uWaveFirstOrder)


thetaOutSecondOrder, thetaOutSecondOrderUncertainty = calcThetaOutAndUncertainty(thetaSecondOrder, THETA_IN, UNCERTAINTY_ANGLE, UNCERTAINTY_IN)
deltaSecondOrder, deltaSecondUncertainty = deltaTotalAndUncertainty(THETA_IN, GRATING_CONSTANT, thetaOutSecondOrder, UNCERTAINTY_IN, thetaOutSecondOrderUncertainty)
wavelengthSecondOrder, uWaveSecondOrder = wavelengthAndUncertainty(2, deltaSecondOrder, deltaSecondUncertainty)
print()
print("Wavelength for Second Order")
print(wavelengthSecondOrder)

print()
print("Uncertainty Wavelength")
print(uWaveSecondOrder)



energyFirstOrder, uEnergyFirst = energyAndUncertainty(wavelengthFirstOrder, uWaveFirstOrder)
print()
print("Total Energy of first order")
print(energyFirstOrder)
print()

print("Uncertainty Total Energyfirst order")
print(uEnergyFirst)

energySecondOrder, uEnergySecond = energyAndUncertainty(wavelengthSecondOrder, uWaveSecondOrder)
print()
print("Total Energy of second order")
print(energySecondOrder)
print()

print("Uncertainty Total Energy second order")
print(uEnergySecond)

print()
print("Wavelength first order")
print(wavelengthFirstOrder)
print()
print("Wavelength Second order")
print(wavelengthSecondOrder)

print()
print("Uncertainty Wavelength first")
print(uWaveFirstOrder)
print()

inverseWaveFirst = 1.0 / wavelengthFirstOrder
inverseWaveFirst = np.array(sorted(list(inverseWaveFirst)))
print(inverseWaveFirst)
#%%
B.pl.clf()
nOne = 1
slopeIn = []
slopeUncertaintyIn = []
slopeOut = []
slopeUncertaintyOut = []

for i in range(3):
    nTwo = np.array([nOne + 1, nOne + 2, nOne + 3, nOne + 4])
    xAxis = (1 / (nOne ** 2)) - (1.0 / (nTwo ** 2))
    
    xAxis = np.array(sorted(list(xAxis)))
    
    uY = inverseLambdaUncertainty(inverseWaveFirst, uWaveFirstOrder)
    
    B.plot_exp(x = xAxis, y = inverseWaveFirst, dy = uY)
    line = B.linefit(x = xAxis, y = inverseWaveFirst)
    slopeIn.append(line.slope)
    slopeUncertaintyIn.append(line.sigma_s)
    
    nOne += 1


B.pl.xlabel(r"$\frac{1}{n_{1}^{2}} - \frac{1}{n_{2}^{2}}$", fontsize = 16)
B.pl.ylabel("Inverse Wavelength (1/m)", fontsize = 16)
B.pl.title("First Order Inverse Wavelength vs " + r"$\frac{1}{n_{1}^{2}} - \frac{1}{n_{2}^{2}}$", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)


inverseWaveSecond = 1.0 / wavelengthSecondOrder
inverseWaveSecond = np.array(sorted(list(inverseWaveSecond)))

nOne = 1
print()

B.pl.clf()
for i in range(3):
    nTwo = np.array([nOne + 1, nOne + 2, nOne + 3, nOne + 4])
    xAxis = (1 / (nOne ** 2)) - (1.0 / (nTwo ** 2))
    
    xAxis = np.array(sorted(list(xAxis)))
    
    uY = inverseLambdaUncertainty(inverseWaveSecond, uWaveSecondOrder)
    
    B.plot_exp(x = xAxis, y = inverseWaveSecond, dy = uY)
    line = B.linefit(x = xAxis, y = inverseWaveSecond)
    slopeOut.append(line.slope)
    slopeUncertaintyOut.append(line.sigma_s)
    
    nOne += 1



#The one used in the lab report was V2
