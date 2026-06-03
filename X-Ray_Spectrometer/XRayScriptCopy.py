# -*- coding: utf-8 -*-

import numpy as np
import LT.box as B
import statistics
from scipy.optimize import curve_fit

# =========================================================
# ===============  Data Loading Functions =================
# =========================================================

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
    """Uncertainty = sqrt(count)"""
    return np.sqrt(count)

# =========================================================
# ===============  Gaussian + Baseline ====================
# =========================================================

def gaussian(x, amplitude, mean, sigma, baseline):
    """Gaussian with baseline correction"""
    return baseline + amplitude * np.exp(-((x - mean) ** 2) / (2 * sigma ** 2))

# =========================================================
# =====================  Main Code  ========================
# =========================================================

# --- Load data ---
degreesCoarse, countCoarse = getDataCoarse('CoarseAnglesMeasurements.txt')
degreesFine, countFine = getDataFine('FineAngleMeasurements.txt')

# Combine both datasets
all_angles = np.append(degreesCoarse, degreesFine)
all_counts = np.append(countCoarse, countFine)

# Sort by angle
sortAngle, sortCount = zip(*sorted(zip(list(all_angles), list(all_counts))))

# Remove duplicate angles
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

# Compute uncertainties
uncertaintyCount = countUncertainty(unique_counts) / 10.0
unique_counts = unique_counts / 10.0  # Normalize

# Estimate background (roughly constant across data)
background = np.median(unique_counts[unique_counts < np.percentile(unique_counts, 20)])
print(f"Estimated background: {background:.2f} counts/s")

# =========================================================
# ====================  Plot Setup  ========================
# =========================================================

B.pl.clf()
B.plot_exp(x=unique_angles, y=unique_counts, dy=uncertaintyCount,
            elinewidth=1.5, capsize=3, alpha=0.7)

B.pl.xlabel("Angle 2θ (°)", fontsize=16)
B.pl.ylabel("Counts per Second", fontsize=16)
B.pl.title("Counts per Second vs Angle", fontsize=16)
B.pl.tick_params(axis='both', labelsize=16)

# =========================================================
# ====================  Peak Fitting  ======================
# =========================================================

# Peak regions (adjust if needed)
peakBeginning = [28.0, 31.0, 59.0, 65.83, 95.0, 109.0]
peakEnding    = [30.0, 33.0, 60.0, 68.0, 96.2, 111.0]
amplitudeGuess = [101.9, 335.8, 25.2, 79.0, 11.8, 34.8]

for a, b, e in zip(amplitudeGuess, peakBeginning, peakEnding):
    peak = (unique_angles >= b) & (unique_angles <= e)
    angles = unique_angles[peak]
    countsPeak = unique_counts[peak]
    
    if len(angles) < 5:
        continue  # skip if too few points

    # Initial guesses: [amplitude, mean, sigma, baseline]
    mean_guess = np.mean(angles)
    sigma_guess = (e - b) / 4.0
    baseline_guess = background

    p0 = [a, mean_guess, sigma_guess, baseline_guess]

    try:
        popt, pcov = curve_fit(gaussian, angles, countsPeak, p0=p0)
        fit_amp, fit_mean, fit_sigma, fit_baseline = popt

        print(f"Fitted peak near {fit_mean:.2f}°:")
        print(f"  Amplitude = {fit_amp:.2f}")
        print(f"  Sigma = {fit_sigma:.4f}")
        print(f"  Baseline = {fit_baseline:.2f}\n")

        # Plot fitted Gaussian
        fine_x = np.linspace(b, e, 300)
        fine_y = gaussian(fine_x, *popt)
        B.plot_line(x=fine_x, y=fine_y, color='r')

    except RuntimeError:
        print(f"Fit failed for region {b}–{e}°")

B.pl.show()
