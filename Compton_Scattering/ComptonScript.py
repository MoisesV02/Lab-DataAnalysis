# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 22:12:27 2025

@author: varga
"""

import numpy as np
import LT.box as B

# get a selection using a plot
def peakpos(h):
   """
   return a selection for data currently plotted

   x array of data that were plotted
   """
   x1,x2 = B.pl.xlim()
   h.fit(x1, x2)
   h.plot_fit()
   return h.mean.value, h.mean.err

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

def electronMass(E0, Ef, theta, uEf):
    masses = []
    uMasses = []
    
    UNCERTAINTY_DEGREE = 1 / 2.0
    UNCERTAINTY_DEGREE = UNCERTAINTY_DEGREE * np.pi / 180.0
    
    
    for a, e, uE in zip(theta, Ef, uEf):
        denom = E0 - e
        numerator = E0 * e * (1 - np.cos(a))
        masses.append(numerator / denom)
        
        
        firstTerm = (1 - np.cos(a)) * ((E0 * (E0 - e)) + (E0 * e)) / ((E0 - e) ** 2)
        secondTerm = E0 * e * np.sin(a) / (E0 - e)
        
        uMasses.append(np.sqrt(((firstTerm * uE) ** 2) + ((secondTerm * UNCERTAINTY_DEGREE) ** 2)))
        
    return np.array(masses), np.array(uMasses)

def uncertaintyGraph(theta, Ef, uEf):
    UNCERTAINTY_DEGREE = 1 / 2.0
    UNCERTAINTY_DEGREE = UNCERTAINTY_DEGREE * np.pi / 180.0
    
    uX = []
    uY = []
    
    for a in theta:
        uX.append(UNCERTAINTY_DEGREE * np.sin(a))
        
    for e, uE in zip(Ef, uEf):
        uY.append(uE / (e ** 2))
        
    return np.array(uX), np.array(uY)
        


#Main

B.pl.clf()
cal = B.get_spectrum('Calibration.Spe')
cal.plot()
B.pl.xlim(0, 550)

#%%
B.pl.clf()
peakPositions = np.array([(np.float64(320.4581680193087), np.float64(0.6635515584879588)),
                 (np.float64(281.84491147471573), np.float64(0.4087387203510716)),
                 (np.float64(163.18290455631254), np.float64(0.18638152291696855))
    ])


pos = np.array(peakPositions[:,0])

energies = np.array([1.3325, 1.1732, 0.6617])       #Energies of known peaks in MeV


B.plot_exp(x = pos, y = energies)
cf = B.linefit(x= pos, y= energies)

B.pl.ylabel('Energies (MeV)', fontsize= 16)
B.pl.xlabel("Positions", fontsize = 16)
B.pl.title("Energies (MeV) vs Positions", fontsize = 16)
B.pl.tick_params(axis='both', labelsize=16)

print()
print("Slope: " + str(cf.slope))
print("Slope Uncertainty: " + str(cf.sigma_s))
print("Chi Squared: " + str(cf.chi_red))
print()

#%%
B.pl.clf()
calNew = B.get_spectrum('Calibration.Spe', cf)
calNew.plot()
B.pl.xlabel("Energies (MeV)")
B.pl.title("Cs-137 and Co-60 Peaks")

B.pl.xlim(0, 2.25)

#%%
B.pl.clf()
sp20_t = B.get_spectrum('20degreeWITHtarget.Spe', cf)
sp20_e = B.get_spectrum('20degreeNOtarget.Spe', cf)


sp20 = sp20_t - sp20_e
sp20 = sp20.rebin(5)

sp20.plot()

B.pl.xlim(0, 2.2)
B.pl.ylim(0, 650)

B.pl.xlabel("Energies (MeV)")
B.pl.ylabel("Counts")
B.pl.title("Photo Peak at 20 Degrees")


#%%
B.pl.clf()
sp40_t = B.get_spectrum('40degreeWITHtarget.Spe', cf)
sp40_e = B.get_spectrum('40degreesNOtarget.Spe', cf)


sp40 = sp40_t - sp40_e
sp40 = sp40.rebin(5)

sp40.plot()

B.pl.xlim(0, 2.2)
B.pl.ylim(0, 500)

B.pl.xlabel("Energies (MeV)")
B.pl.ylabel("Counts")
B.pl.title("Photo Peak at 40 Degrees")

#%%
B.pl.clf()
sp60_t = B.get_spectrum('60degreeWITHtarget.Spe', cf)
sp60_e = B.get_spectrum('60degreeNOtarget.Spe', cf)


sp60 = sp60_t - sp60_e
sp60 = sp60.rebin(5)

sp60.plot()

B.pl.xlim(0, 2.2)
B.pl.ylim(0, 850)

B.pl.xlabel("Energies (MeV)")
B.pl.ylabel("Counts")
B.pl.title("Photo Peak at 60 Degrees")

#%%
B.pl.clf()
sp80_t = B.get_spectrum('80degreeWITHtarget.Spe', cf)
sp80_e = B.get_spectrum('80degreeNOtarget.Spe', cf)


sp80 = sp80_t - sp80_e
sp80 = sp80.rebin(4)

sp80.plot()

B.pl.xlim(0, 2.2)
B.pl.ylim(0, 800)

B.pl.xlabel("Energies (MeV)")
B.pl.ylabel("Counts")
B.pl.title("Photo Peak at 80 Degrees")

#%%
B.pl.clf()
sp100_t = B.get_spectrum('100degreesWITHtarget.Spe', cf)
sp100_e = B.get_spectrum('100degreeNOtarget.Spe', cf)


sp100 = sp100_t - sp100_e
sp100 = sp100.rebin(4)

sp100.plot()

B.pl.xlim(0, 2.2)
B.pl.ylim(0, 900)

B.pl.xlabel("Energies (MeV)")
B.pl.ylabel("Counts")
B.pl.title("Photo Peak at 100 Degrees")

#%%
B.pl.clf()
sp120_t = B.get_spectrum('120degreeWITHtarget.Spe', cf)
sp120_e = B.get_spectrum('120degreeNOtarget.Spe', cf)


sp120 = sp120_t - sp120_e
sp120 = sp120.rebin(2)

sp120.plot()

B.pl.xlim(0, 2.2)
B.pl.ylim(0, 800)

B.pl.xlabel("Energies (MeV)")
B.pl.ylabel("Counts")
B.pl.title("Photo Peak at 120 Degrees")

#%%
energyPeaks = np.array([(np.float64(0.6181337696229872), np.float64(0.05720830847112319)),
                    (np.float64(0.48066392702302246), np.float64(0.00934025128133461)),
                    (np.float64(0.39244515431759286), np.float64(0.008066495690872314)),
                    (np.float64(0.30594821430753355), np.float64(0.003016320001956094)),
                    (np.float64(0.24697340633666542), np.float64(0.0031164487466001285)),
                    (np.float64(0.21196216256704214), np.float64(0.0012010248137051393))
                    ])

B.pl.clf()


energyFinal = energyPeaks[:,0]
uEnergyFinal = energyPeaks[:,1]

theta = np.array([20, 40, 60, 80, 100, 120])
thetaRad = theta * np.pi / (180.0)

xAxis = 1 - np.cos(thetaRad)

uX, uY = uncertaintyGraph(thetaRad, energyFinal, uEnergyFinal)

B.plot_exp(x = xAxis, y = (1.0 / energyFinal), dy = uY, xerr = uX)
line = B.linefit(x = xAxis, y= (1.0 / energyFinal), yerr = uY)

B.pl.xlabel(r"$1 - \cos\theta$", fontsize=16)
B.pl.ylabel(r"$\frac{1}{E_f}$ (MeV$^{-1}$)", fontsize=16)
B.pl.title(r"$\frac{1}{E_f} \ \mathrm{vs.} \ 1 - \cos\theta$", fontsize=16)
B.pl.tick_params(axis='both', labelsize=16)


slope = line.slope
uSlope = line.sigma_s

print()
print("Slope: " + str(slope))
print("Uncertainty Slope: " + str(uSlope))
print("Chi Square: " + str(line.chi_red))

mass_e = 1.0 / slope
uMass_e = uSlope / (slope ** 2)

print()
print(f"Experimental slope mass of electron (MeV): {mass_e:.3} ± {uMass_e:.3}")



electronMass, uElectronMass = electronMass(energies[2], energyFinal, thetaRad, uEnergyFinal)

avgMass, avgUncertainty = calcMeanandUncertainty(electronMass, uElectronMass)

print()
print(f"Experimental mass of electron (MeV): {avgMass:.3} ± {avgUncertainty:.3}")