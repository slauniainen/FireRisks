# -*- coding: utf-8 -*-
"""
Runs few test scenarios and draws figures

Created on Mon Jun 21 09:35:42 2021

@author: 03081268
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit

from moss_cyano import fDepo, fT, fWater
from mcn2fix import driver
from utils import read_sitedata, read_forcing

EPS = np.finfo(float).eps

# read forcing
rcp_scen = 'rcp45'
#grid_id = 3812 #[68.6N / 25.2E]
#grid_id = 3787 #[65.6N / 25.6E]
#grid_id = 3286 #[63.5N / 24.4E]
grid_id = 3269 #[61.8N / 24.4E]

forc_present = read_forcing(rcp_scen, grid_id, fdate='1990-01-01', ldate='2020-12-31')
forc_future = read_forcing(rcp_scen, grid_id, fdate='2069-01-01', ldate='2099-12-31') 
# parameters

# canopy waterbudget
water_p = {'interception': {'wmax': 0.2, 'alpha': 1.28},
          'snowpack': {'Tm': 0.0, 'Km': 2.9e-5, 'Kf': 5.8e-6, 'R': 0.05, 
                       'Tmax': 0.5, 'Tmin': -0.5, 'Sliq': 0.0, 'Sice': 0.0},
          'organic_layer': {'Wmax': 10.0, 'Wmin': 0.1, 'Wcrit': 4.0,
                            'alpha': 1.28, 'Wliq': 10.0}                  
          }

# read sitedata sfile = r'Data/ICP_L2_sitedata.csv'
          
sfile = r'Data/ICP_L2_sitedata.csv'
#siteid = 3 # Pallas P
sideid = 4 # Pallas S
sitedata = read_sitedata(sfile, siteid=4)

nitrogen_p = { 'Hylo': {'potrate': [25.6, 12.8], 'DM': [sitedata['Hylo_u'], sitedata['Hylo_l']], 'wresponse': [0.44, 7.0]},
               'Pleu': {'potrate':[27.0, 7.3], 'DM': [sitedata['Pleu_u'] + sitedata['Other_u'], sitedata['Pleu_l'] + sitedata['Other_l']], 
                        'wresponse': [0.44, 7.0]},
               'Dicr': {'potrate': [0.2, 0.0], 'DM': [sitedata['Dicr_u'], sitedata['Dicr_l']], 'wresponse': [0.44, 7.0]},
             }

#%% Run model for pallas P

#'date', 'N2fix', 'cN2fix', 'water_content', 'swe', 'fT', 'fw'
               
res0, mod0 = driver(sitedata, forc_present, nitrogen_p, water_p)
res1, mod1 = driver(sitedata, forc_future, nitrogen_p, water_p)

doy0 = forc_present.doy.values
doy1 = forc_future.doy.values

#%% group data to annual balances

res0['fclim'] = res0['fT']*res0['fw']
res1['fclim'] = res1['fT']*res1['fw']

r0 = res0.groupby('doy')
r1 = res1.groupby('doy')
f0 = forc_present.groupby('doy')
f1 = forc_future.groupby('doy')

doy0 = r0.mean().index.values
doy1 = r1.mean().index.values

fig, ax = plt.subplots(4,1)
fig.set_size_inches(11.5, 11.5)

ax[0].fill_between(doy0, 1e3*r0.quantile(0.25)['N2fix'], 1e3*r0.quantile(0.75)['N2fix'], color='b', alpha=0.3)
ax[0].plot(doy0, 1e3*r0.mean()['N2fix'], 'b-', label='1990-2020')

ax[0].fill_between(doy1, 1e3*r1.quantile(0.25)['N2fix'], 1e3*r1.quantile(0.75)['N2fix'], color='r', alpha=0.3)
ax[0].plot(doy1, 1e3*r1.mean()['N2fix'], 'r-', label='2070-2100')
ax[0].legend()
ax[0].set_ylabel('N$_F$ (g N ha-1 d-1)'); ax[0].set_xlabel('doy')

ax[1].fill_between(doy0, r0.quantile(0.25)['fT'], r0.quantile(0.75)['fT'], color='b', alpha=0.3)
ax[1].plot(doy0, r0.mean()['fT'], 'b-', label='1990-2020')

ax[1].fill_between(doy1, r1.quantile(0.25)['fT'], r1.quantile(0.75)['fT'], color='r', alpha=0.3)
ax[1].plot(doy1, r1.mean()['fT'], 'r-', label='2070-2100')
ax[1].set_ylabel('fT (-)')

ax[2].fill_between(doy0, r0.quantile(0.25)['fw'], r0.quantile(0.75)['fw'], color='b', alpha=0.3)
ax[2].plot(doy0, r0.mean()['fw'], 'b-', label='1990-2020')

ax[2].fill_between(doy1, r1.quantile(0.25)['fw'], r1.quantile(0.75)['fw'], color='r', alpha=0.3)
ax[2].plot(doy1, r1.mean()['fw'], 'r-', label='2070-2100')
ax[2].set_ylabel('fw (-)')

ax[3].fill_between(doy0, r0.quantile(0.25)['fclim'], r0.quantile(0.75)['fclim'], color='b', alpha=0.3)
ax[3].plot(doy0, r0.mean()['fclim'], 'b-', label='1990-2020')

ax[3].fill_between(doy1, r1.quantile(0.25)['fclim'], r1.quantile(0.75)['fclim'], color='r', alpha=0.3)
ax[3].plot(doy1, r1.mean()['fclim'], 'r-', label='2070-2100')
ax[3].set_ylabel('fclim (-)')
ax[3].set_xlabel('doy')

plt.savefig(r'Data/Pallas_pine_Nfix_IAV.png', dpi=900)

fig, ax = plt.subplots(4,1)
fig.set_size_inches(11.5, 11.5)

ax[0].fill_between(doy0, f0.quantile(0.25)['Tair'], f0.quantile(0.75)['Tair'], color='b', alpha=0.3)
ax[0].plot(doy0, f0.mean()['Tair'], 'b-', label='1990-2020')
ax[0].set_ylabel('Ta (degC)')
        
ax[0].fill_between(doy1, f1.quantile(0.25)['Tair'], f1.quantile(0.75)['Tair'], color='r', alpha=0.3)
ax[0].plot(doy1, f1.mean()['Tair'], 'r-', label='2070-2100')
ax[0].legend()

ax[1].fill_between(doy0, f0.quantile(0.25)['Prec']*86400, f0.quantile(0.75)['Prec']*86400, color='b', alpha=0.3)
ax[1].plot(doy0, f0.mean()['Prec']*86400, 'b-', label='1990-2020')
ax[1].set_ylabel('Prec (mm d-1)')

ax[1].fill_between(doy1, f1.quantile(0.25)['Prec']*86400, f1.quantile(0.75)['Prec']*86400, color='r', alpha=0.3)
ax[1].plot(doy1, f1.mean()['Prec']*86400, 'r-', label='2070-2100')

ax[2].fill_between(doy0, r0.quantile(0.25)['swe'], r0.quantile(0.75)['swe'], color='b', alpha=0.3)
ax[2].plot(doy0, r0.mean()['swe'], 'b-', label='1990-2020')
ax[2].set_ylabel('SWE (mm)')
ax[2].fill_between(doy1, r1.quantile(0.25)['swe'], r1.quantile(0.75)['swe'], color='r', alpha=0.3)
ax[2].plot(doy1, r1.mean()['swe'], 'r-', label='2070-2100')

#plt.savefig(r'Data/Pallas_pine_climate.png', dpi=900))

            #%% annual histograms

a0 = res0.resample('A')
a1 = res1.resample('A')

fig, ax = plt.subplots(1)

c='b'
ax.boxplot(a0['N2fix'].sum(), 'r', sym='o', positions=[1], 
           boxprops=dict(color=c),
           capprops=dict(color=c),
           whiskerprops=dict(color=c),
           flierprops=dict(color=c, markeredgecolor=c),
           medianprops=dict(color=c)
           )
c='r'
ax.boxplot(a1['N2fix'].sum(), 'b', sym='o', positions=[1.5],
           boxprops=dict( color=c),
           capprops=dict(color=c),
           whiskerprops=dict(color=c),
           flierprops=dict(color=c, markeredgecolor=c),
           medianprops=dict(color=c)
           )
ax.set_ylabel('F$_N$ (kg ha-1 a-1)')
ax.set_xticklabels(['current', 'future'])


capa = 1.85 / 1000 / 1000 #kgNkgdm-1

DM = sum(nitrogen_p['Hylo']['DM'] + nitrogen_p['Pleu']['DM'] + nitrogen_p['Hylo']['DM'])
NS = capa * DM  * 200

plt.plot([0.5, 2], [NS, NS], 'k--')
#plt.savefig(plt.savefig(r'Data/N2fix_Pallas_pine_current_future.png', dpi=900)

#fig, ax = plt.subplots(1,2)
#
#d0 = np.histogram(a0['N2fix'].sum().values)
#x0 = d0[1][1:]
#dm0 = np.mean(a0['N2fix'].sum().values)
#d1 = np.histogram(a1['N2fix'].sum().values)
#x1 = d1[1][1:]
#dm1 = np.mean(a1['N2fix'].sum().values)
#
#ax[0].bar(x0, d0[0]/max(d0[0]), color='b', alpha=0.5, width=d0[1][0]/10, label='1990-2020')
#ax[0].plot([dm0, dm0], [0,1], 'bo-')
#ax[0].bar(x1, d1[0]/max(d1[0]), color='r', alpha=0.5, width=d1[1][0]/10, label='2070-2100')
#ax[0].plot([dm1, dm1], [0,1], 'ro-')

