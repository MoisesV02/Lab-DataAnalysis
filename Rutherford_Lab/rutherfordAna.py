#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 12:00:46 2025

@author: hszumila
"""

import numpy as np
import LT.box as B
import matplotlib.pyplot as plt

def get_time(sp):
       tt = float(sp.title.split('live time:')[1].split('s,')[0].strip())
       return tt
   
def S(x):                       # define the function
       return C()/B.np.sin(0.5*(x - t0()))**4

#read the input files
spec0=B.get_spectrum('0_deg.Spe')
spec0e=B.get_spectrum('0_deg_empty.Spe')
spcfiles=["5_deg.Spe","n5_deg.Spe","10_deg.Spe","n10_deg.Spe","15_deg.Spe",
          "n15_deg.Spe","20_deg.Spe","n20_deg.Spe","n30_deg.Spe","n40_deg.Spe"]
angles=[5,-5,10,-10,15,-15,20,-20,-30,-40]

#plot the 0 deg empty and data on a single plot with two different colors
#%%
spec0.plot()
spec0e.plot()
B.pl.show()
B.pl.title("Spectrum With & Without Gold Sheet at 0 Degrees")
#report the statistical uncertainty, count rates with each angle
#%%
B.pl.clf()
rate=np.array([])
rateU=np.array([])
theta_r=np.array([])
for ii in range(len(spcfiles)):
    #B.pl.clf()
    data=B.get_spectrum(spcfiles[ii])
    sp=data.plot()
    C, dC = data.sum(50, 100) # get the events
    t = get_time(data)   # get the time
    R = C/t            # calculate the rate
    dR = dC/t          # calculate the error in the rate
    rate=np.append(rate,R)
    rateU=np.append(rateU,dR)
    theta_r=np.append(theta_r,angles[ii]*(22./7.)/180.0)
    print('Angle:',angles[ii],"deg,has ",C," events, N/time: ",f"{R:.2f}","+/-",f"{dR:.2f}","over ",t," s.")
#%%
#plot the ln(count rate) vs sin(theta/2) with uncertainty, abs val of theta
B.pl.clf()
B.plot_exp(theta_r, rate, rateU, logy=True,x_label = 'theta [rad]', 
           y_label = 'ln(count rate)', plot_title = 'Log count rate vs angle')
#%%

#determine an offset in the angle to get the points aligned
C = B.Parameter(10., 'C')       # define the parameter amplitude  
t0 = B.Parameter(0, 'theta_0') # define parameter for offset
sfit = B.genfit(S, [C, t0], x = theta_r, y = rate,y_err = rateU) 
B.pl.show()

offset = -0.004068
aTheta = abs(theta_r-offset)/2.
#lth = aTheta >0.0 
xvar = np.log(np.sin(aTheta))
mAng = xvar>-2.3
#%%
#fit the data with a straight line and compare to rutherford
B.pl.clf()
B.plot_exp(xvar, rate, rateU, logy=True, x_label = 'ln(sin(θ/2)', 
           y_label = 'ln(count rate)', plot_title = 'Log count rate vs log angle')
coeffs = np.polyfit(xvar[mAng], np.log(rate[mAng]), 1, w=1/np.log(rateU[mAng]))
a, b = coeffs
print(f"Slope = {a:.3f}, Intercept = {b:.3f}")
R_fit = np.exp(a * xvar + b)
plt.semilogy(xvar, R_fit, '-', label='Fit')
B.pl.show()


print()
print("C =", C(), "+/-", C.err)
print("theta_0 =", t0(), "+/-", t0.err)


