# -*- coding: utf-8 -*-
"""
Sandbox for Moss - cyanobacteria Nfix
Created on Thu Jun 17 13:13:12 2021

@author: 03081268
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from moss_cyano import fDepo, fT, fWater
EPS = np.finfo(float).eps


#%% massage Salemaa et al. 2019 data.
    
dat = pd.read_csv(r'Data\Salemaa_etal_2019Table2.txt', sep=';', usecols=[1,2,3,4])

#%%
#pot = pd.read_csv(r'Data\Salemaa_etal_2019TableA2.txt', sep=';', skiprows=0)
#hylo = pot.iloc[0:5]
#pleu = pot.iloc[5:]

#plt.figure()
#
#for k in ['Dicr', 'Hylo', 'Pleu']:
#    plt.plot(dat['Ndepo'], dat[k].values / np.nanmax(dat[k].values), 'o', label=k)
#plt.legend()

N = 3 * list(dat['Ndepo'].values)
C = []
for k in ['Dicr', 'Hylo', 'Pleu']:
    C.extend(dat[k].values / np.nanmax(dat[k].values))

D = pd.DataFrame(columns=['N', 'f'], data=None)
D['N'] = N; D['f'] = C
D = D.dropna(axis=0)

# global fit
popt0, pcov0 = curve_fit(fDepo, D['N'].values, D['f'].values, bounds=(-5., 0.))

perr0 = np.sqrt(np.diag(pcov0))

xx = np.arange(0.5, 4.5, 0.1)


plt.figure()
c = plt.rcParams['axes.prop_cycle'].by_key()['color']
n = 0

for k in ['Dicr', 'Hylo', 'Pleu']:
    plt.plot(dat['Ndepo'], dat[k] / np.nanmax(dat[k].values), 'o', color=c[n])
    x = dat['Ndepo'].values
    y = dat[k].values / np.nanmax(dat[k].values)
    ix = np.where(y>=0)
    x = x[ix]
    y = y[ix]
    
    # fit per species
    popt, pcov = curve_fit(fDepo, x, y, bounds=(-5., 0.))
    fit = fDepo(xx, popt[0])
    plt.plot(xx, fit, '--', color=c[n], label=k + ' exp[%.2f (x-1)]' % (popt[0]))
    n += 1

fit = fDepo(xx, popt0[0])
plt.plot(xx, fit, 'k-', label='exp[%.2f (x-1)]' % (popt0[0]))
plt.legend(fontsize=10)
plt.xlabel('N deposition (kg ha$^{-1}$ a$^{-1}$)')
plt.ylabel('fN$_{d}$ (-)')
plt.savefig(r'Data/fPot.png', dpi=900)

#%% Plot Temperature and water modifiers

t = np.arange(0, 35, 0.1)
w = np.arange(0, 10, 0.1)
ndepo = np.arange(0.5, 4.0, 0.1)

plt.figure()
plt.subplot(221); plt.plot(ndepo, fDepo(ndepo), 'k-'); plt.plot(1.0, fDepo(1.0), 'ko'); plt.ylabel('fN$_d$ (-)'); plt.xlabel('N$_d$ (kg ha-1 a-1)')
plt.subplot(222); plt.plot(t, fT(t), 'k-'); plt.plot(25.0, fT(25.0), 'ko'); plt.ylabel('fT (-)'); plt.xlabel('T (degC)')
plt.subplot(223); plt.plot(w, fWater(w), 'k-'); plt.plot(10.0, 1.0, 'ko'); plt.ylabel('fW (-)'); plt.xlabel('W (g g-1)')
plt.savefig(r'Data/fNd_fT_fW.png', dpi=900)

#%%