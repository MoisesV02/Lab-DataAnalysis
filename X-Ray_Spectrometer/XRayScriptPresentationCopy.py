# -*- coding: utf-8 -*-
#Used this version for the lab report
"""
Created on Fri Sep 26 11:09:37 2025

@author: varga
"""

import numpy as np
import LT.box as B
import statistics

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


def getDataCoarse(fileName):
    fileObj = B.get_file(fileName)
    
    degrees = fileObj['Angle']
    count = fileObj['Count']
    
    return degrees, count

def getDataFine(fileName):
    fileObj = B.get_file(fileName)
    
    degrees = fileObj['Angle']
    count = fileObj['Count']
    arcMin = fileObj['ArcMinute']
        
    degrees = degrees + (arcMin / 60.0)
    
    return degrees, count

def countUncertainty(count):
    uncertainty = []
    
    for c in count:
       uncertainty.append(np.sqrt(c)) 
    
    
    return np.array(uncertainty)

def gaussian(array, amplitude):
    probabilities = []
    mean = statistics.mean(array)
    variance = statistics.variance(array)
    for value in array:
        exponent = - ((value - mean) ** 2) / (2 * variance)
        product = amplitude * np.exp(exponent)
        probabilities.append(product)
     
    print("Mean: " + str(mean))
    return np.array(probabilities)

def uncertaintySin(angles):
    uAngle = 0.58 * np.pi / 180.0
    uncertainty = []
    for a in angles:
        uncertainty.append( np.abs(uAngle * np.cos(a) / 2.0))
        
    return np.array(uncertainty)


#Main
degreesCoarse, countCoarse = getDataCoarse('CoarseAnglesMeasurements.txt')
degreesFine, countFine = getDataFine('FineAngleMeasurements.txt')

all_angles = np.append(degreesCoarse, degreesFine)


all_counts = np.append(countCoarse, countFine)


sortAngle, sortCount = zip(*sorted(zip(list(all_angles), list(all_counts))))


unique_angles = []
unique_counts = []
seen_angles = set()

for i in range(len(sortAngle)):
    if sortAngle[i] not in seen_angles:
        unique_angles.append(sortAngle[i])
        unique_counts.append(sortCount[i])
        seen_angles.add(sortAngle[i])
        
unique_angles = np.array(unique_angles)
unique_counts = np.array(unique_counts)
uncertaintyAngle = np.array([0.58] * len(unique_angles))

uncertaintyCount = countUncertainty(unique_counts) / 10.0

unique_counts = unique_counts / 10.0


B.pl.clf()


B.plot_exp(x = unique_angles, y = unique_counts, dy = uncertaintyCount, xerr = uncertaintyAngle, elinewidth= 1.5, capsize= 3, alpha = 0.7)


B.pl.xlabel("Angle 2θ (°)", fontsize = 18)
B.pl.ylabel("Counts per Second", fontsize = 18)
B.pl.title("Counts per Second vs Angle", fontsize = 18)
B.pl.tick_params(axis='both', labelsize=18)


peakBeginning = [28.0, 31.0, 59.0, 65.83, 95.0, 109.0]
peakEnding = [30.0, 33.0, 60.0, 68.0, 96.2, 111.0]
amplitude = [101.9, 335.8, 25.2, 79.0, 11.8, 34.8]


for a, b, e in zip(amplitude, peakBeginning, peakEnding):
    peak = np.array(unique_angles >= b) & (unique_angles <= e)
    angles = unique_angles[peak]
    countsPeak = unique_counts[peak]
    
    probabilities = gaussian(angles, a)
    
    B.plot_line(x = angles, y = probabilities, color='r')
    
B.pl.annotate(r"$K_{\beta_1}$", (22, 96), color='b', fontsize = 16)
B.pl.annotate(r"${29.18^\circ}$", (21, 75), color = 'b', fontsize = 12)

B.pl.annotate(r"$K_{\alpha_1}$", (34, 333), color='b', fontsize = 16)
B.pl.annotate(r"${32^\circ}$", (34.3, 313), color = 'b', fontsize = 12)



B.pl.annotate(r"$K_{\beta_2}$", (55, 47), color='b', fontsize = 16)
B.pl.annotate(r"${59.5^\circ}$", (54, 27), color = 'b', fontsize = 12)

B.pl.annotate(r"$K_{\alpha_2}$", (69, 92), color='b', fontsize = 16)
B.pl.annotate(r"${66.5^\circ}$", (68.5, 72), color = 'b', fontsize = 12)



B.pl.annotate(r"$K_{\beta_3}$", (89, 36), color='b', fontsize = 16)
B.pl.annotate(r"${95.3^\circ}$", (88, 16), color = 'b', fontsize = 12)

B.pl.annotate(r"$K_{\alpha_3}$", (113.5, 45), color='b', fontsize = 16)
B.pl.annotate(r"${110.2^\circ}$", (112.5, 25), color = 'b', fontsize = 12)

#%%
B.pl.clf()
peakAngles = [unique_angles[11], unique_angles[18], unique_angles[50], unique_angles[63], unique_angles[96], unique_angles[117]]
peakAngles = np.array(peakAngles)

peakAngles = peakAngles * np.pi / 180.0
peakAngles = peakAngles / 2.0

#Beta
peakAnglesFirstOrder = np.array([peakAngles[0], peakAngles[2], peakAngles[4]])

#Alpha
peakAnglesSecondOrder = np.array([peakAngles[1], peakAngles[3], peakAngles[5]])


n = np.array([1,2,3])

uY = uncertaintySin(peakAngles)
uYFirstOrder = np.array([uY[0], uY[2], uY[4]])
uYSecondOrder = np.array([uY[1], uY[3], uY[5]])

B.plot_exp(x = n, y = np.sin(peakAnglesFirstOrder), dy = uYFirstOrder)
lineFirstOrder = B.linefit(x = n, y = np.sin(peakAnglesFirstOrder), yerr = uYFirstOrder)

B.plot_exp(x = n, y = np.sin(peakAnglesSecondOrder), dy = uYSecondOrder)
lineSecondOrder = B.linefit(x = n, y = np.sin(peakAnglesSecondOrder), yerr = uYSecondOrder)

B.pl.xlabel("Order of Diffraction (n)", fontsize = 16)
B.pl.ylabel("sin(θ/2) in radians", fontsize = 16)
B.pl.title("sin(θ/2) vs Order of Diffraction", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)

B.pl.annotate(r"$K_{\beta}$", (2.3, 0.53), color='orange', fontsize = 18)
B.pl.annotate(r"$K_{\alpha}$", (2.5, 0.72), color='r', fontsize = 18)

#Alpha line > beta line, alpha is red while beta orange

print()

slopeFirstOrder = lineFirstOrder.slope
uSlopeFirstOrder = lineFirstOrder.sigma_s

slopeSecondOrder = lineSecondOrder.slope
uSlopeSecondOrder = lineSecondOrder.sigma_s
#wavelength = [1.542E-10, 1.392E-10]
wavelength = [1.392E-10, 1.542E-10]

distance = []
uDistance = []
slopes = [slopeFirstOrder, slopeSecondOrder]
uSlopes = [uSlopeFirstOrder, uSlopeSecondOrder]
for w, slope, uSlope in zip(wavelength, slopes, uSlopes):
    d = w / (2 * slope)
    print()
    print(d * 1E9)
    distance.append(d * 1E9)
    
    print((w  * uSlope / (2 * (slope ** 2))) * 1E9)
    uDistance.append((w  * uSlope / (2 * (slope ** 2))) * 1E9)

print()
    
distanceAvg, uAvg = calcMeanandUncertainty(distance, uDistance)

print(f"Experimental Interplanar Spacing: {distanceAvg:.3} ± {uAvg:.3}")