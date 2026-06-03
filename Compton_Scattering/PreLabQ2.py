# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 22:15:03 2025

@author: varga
"""

import numpy as np
#Determine the photon energy if the photon is scattered at 30, 40, 60, and 80 degrees. Do not use calculators
#or any sort. Take a screen shot of how you did it in the python console.

E_o = 0.6617 #MeV
restEnergyElectron = 0.511 #MeV

angles = np.array([30, 40, 60, 80])

angles = angles * np.pi / 180.0

for a in angles:
    photonEnergy = E_o / (1  + E_o / restEnergyElectron * (1 - np.cos(a)))
    
    print(f"At {a * 180.0 / np.pi:2.0f} degrees, energy of photon is: {photonEnergy:.3f} MeV")
