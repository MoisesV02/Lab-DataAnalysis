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
    sigma = np.sqrt(variance)  # Standard deviation of Gaussian
    for value in array:
        exponent = - ((value - mean) ** 2) / (2 * variance)
        product = amplitude * np.exp(exponent)
        probabilities.append(product)
    print("Mean: " + str(mean))
    print("Sigma: " + str(sigma))
    return np.array(probabilities), mean, sigma  # Return mean and sigma

def uncertaintySin(angles, uPeak):
    uAngle = 0.17 * np.pi / 180.0  # Change angle uncertainty to 0.17° and convert to radians
    uncertainty = []
    for a, u in zip(angles, uPeak):
        uncertainty.append(np.abs(np.sqrt((uAngle ** 2) + (u ** 2)) * np.cos(a) / 2.0))
    return np.array(uncertainty)


# Main
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
uncertaintyAngle = np.array([0.17] * len(unique_angles))  # Angle uncertainty is now 0.17°
uncertaintyCount = countUncertainty(unique_counts) / 10.0
unique_counts = unique_counts / 10.0


B.pl.clf()

# Plotting the raw data
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
    
    probabilities, mean, sigma = gaussian(angles, a)  # Get sigma from Gaussian fit
    # Print the Gaussian fit parameters
    print(f"Gaussian Fit - Mean: {mean}, Sigma: {sigma}")

    # Propagate the uncertainty including Gaussian sigma
    # Assuming the Gaussian uncertainty is propagated similarly to angle uncertainty
    
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

# Peak angles (in radians)
peakAngles = [unique_angles[11], unique_angles[18], unique_angles[50], unique_angles[63], unique_angles[96], unique_angles[117]]
peakAngles = np.array(peakAngles)

peakAngles = peakAngles * np.pi / 180.0  # Convert to radians
peakAngles = peakAngles / 2.0  # Apply Bragg's law (θ/2)

# First and second order peaks
peakAnglesFirstOrder = np.array([peakAngles[0], peakAngles[2], peakAngles[4]])
peakAnglesSecondOrder = np.array([peakAngles[1], peakAngles[3], peakAngles[5]])

uPeaks = np.array([0.012, 0.005, 0.014, 0.008, 0.042, 0.029])
uPeaks = uPeaks * np.pi / 180.0  # Convert to radians

n = np.array([1,2,3])

uY = uncertaintySin(peakAngles, uPeaks)
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

# Calculate slopes and uncertainties
slopeFirstOrder = lineFirstOrder.slope
uSlopeFirstOrder = lineFirstOrder.sigma_s

slopeSecondOrder = lineSecondOrder.slope
uSlopeSecondOrder = lineSecondOrder.sigma_s

# Wavelengths
wavelength = [1.392E-10, 1.542E-10]

distance = []
uDistance = []

slopes = [slopeFirstOrder, slopeSecondOrder]
uSlopes = [uSlopeFirstOrder, uSlopeSecondOrder]

for w, slope, uSlope in zip(wavelength, slopes, uSlopes):
    d = w / (2 * slope)
    distance.append(d * 1E9)  # Convert to nm
    uDistance.append((w * uSlope / (2 * (slope ** 2))) * 1E9)  # Error propagation

# Calculate average distance and uncertainty
distanceAvg, uAvg = calcMeanandUncertainty(distance, uDistance)

print(f"Experimental Interplanar Spacing: {distanceAvg:.3} ± {uAvg:.3}")
