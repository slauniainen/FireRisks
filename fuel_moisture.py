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
from utils import read_forcing

from waterbalance import Canopywaterbudget
EPS = np.finfo(float).eps

def driver(dt, forc, LAI, params):
    
    task = Canopywaterbudget(dt, LAI, params)
    
    N = len(forc)
    
    results = {'date': forc['date'], 'doy': forc['doy'], 
               'Prec': np.ones((N,1))*np.NaN, 
               'throughfall': np.ones((N,1))*np.NaN,
               'swe': np.ones((N,1))*np.NaN,
               'canopy_water_storage': np.ones((N,1))*np.NaN, 
               'litter_water_content': np.ones((N,1))*np.NaN,
               'soil_water_content': np.ones((N,1))*np.NaN,
               'soil_water_potential': np.ones((N,1))*np.NaN,
               'transpiration_rate': np.ones((N,1))*np.NaN,
               'et_rate': np.ones((N,1))*np.NaN,
               'dummy': np.ones((N,1))*np.NaN,
               #'Gc0': np.ones((N,1))*np.NaN,
               }
    
    for k in range(N):
        #print(k/N*100)
        
        Trfall, Ec, Ebl, Tr, Wc, Wbl, Ws, Psi_s = task.run(forc.iloc[k].to_dict())
        
        results['throughfall'][k] = Trfall * task.dt
        results['swe'][k] = task.snow.SWE
        results['et_rate'][k] = (Tr + Ec + Wc) * task.dt
        results['transpiration_rate'][k] = Tr * task.dt
        results['canopy_water_storage'][k] = Wc
        results['litter_water_content'][k] = Wbl
        results['soil_water_content'][k] = Ws
        results['soil_water_potential'][k] = Psi_s
        results['dummy'][k] = task.canopy.fPheno
        #results['Gc0'][k] = Gc0
        
    #for k in results.keys():
    #    results[k] = results[k].ravel()
    #results = pd.DataFrame.from_dict(results)
    #results.index = forc.index
    #results = results.drop(columns={'date'})
    results['doy'] = forc.doy
    return results, task



#%%
# read forcing
rcp_scen = 'rcp45'
#grid_id = 3812 #[68.6N / 25.2E]
#grid_id = 3787 #[65.6N / 25.6E]
#grid_id = 3286 #[63.5N / 24.4E]
grid_id = 3269 #[61.8N / 24.4E]

forc_present = read_forcing(rcp_scen, grid_id, fdate='2000-01-01', ldate='2000-12-31')
#forc_future = read_forcing(rcp_scen, grid_id, fdate='2069-01-01', ldate='2099-12-31') 
# parameters

# canopy waterbudget
params = {'interception': {'LAI': None, 'wmax': 0.2, 'alpha': 1.28},
          'snowpack': {'Tm': 0.0, 'Km': 2.9e-5, 'Kf': 5.8e-6, 'R': 0.05, 
                       'Tmax': 0.5, 'Tmin': -0.5, 'Sliq': 0.0, 'Sice': 0.0},
          'organic_layer': {'DM': 0.25, 'Wmax': 10.0, 'Wmin': 0.1, 'Wcrit': 4.0,
                            'alpha': 1.28, 'Wliq': 10.0},
          'soil': {'depth': 0.4, 'Ksat': 1e-5, 'MaxPond': 0.0,
                   'pF': {'ThetaS': 0.54, 'ThetaR': 0.00, 'alpha': 0.448, 'n':1.2} # mesic sites
                  },
          'transpiration': {'Amax': 9.0, 'g1': 2.0, 'q50': 50.0, 'kp': 0.6,
                            'rw': 0.20, 'rwmin': 0.02,
                            'smax': 18.5, 'tau': 13.0, 'xo': -4.0, 'fmin': 0.05
                            }
          
          }

LAI = 3.0
dt = 86400.0

#%% Run model for test location
               
res0, mod0 = driver(dt, forc_present, LAI, params)


# #%% group data
# res0['crown_sat'] = res0['canopy_water_storage'] / mod0.waterbalance.interception.Wmax
# res0['litter_sat'] = res0['litter_water_content'] / mod0.waterbalance.bottomlayer.Wmax

# res1['crown_sat'] = res1['canopy_water_storage'] / mod1.waterbalance.interception.Wmax
# res1['litter_sat'] = res1['litter_water_content'] / mod1.waterbalance.bottomlayer.Wmax

# res2['crown_sat'] = res2['canopy_water_storage'] / mod2.waterbalance.interception.Wmax
# res2['litter_sat'] = res2['litter_water_content'] / mod2.waterbalance.bottomlayer.Wmax

# r0 = res0.groupby('doy')
# r1 = res1.groupby('doy')
# r2 = res2.groupby('doy')
# f0 = forc_present.groupby('doy')
# f1 = forc_future.groupby('doy')

# doy0 = r0.mean().index.values
# doy1 = r1.mean().index.values

# fig, ax = plt.subplots(3,1)
# fig.set_size_inches(11.5, 11.5)

# ax[0].fill_between(doy0, r0.quantile(0.25)['litter_sat'], r0.quantile(0.75)['litter_sat'], color='b', alpha=0.3)
# ax[0].plot(doy0, r0.median()['litter_sat'], 'b-', label='1990-2020, sparse')

# ax[0].fill_between(doy1, r1.quantile(0.25)['litter_sat'], r1.quantile(0.75)['litter_sat'], color='r', alpha=0.3)
# ax[0].plot(doy1, r1.median()['litter_sat'], 'r-', label='1990-2020, dense')

# ax[0].legend()
# ax[0].set_ylabel('litter saturation'); ax[0].set_xlabel('doy')

# ax[1].fill_between(doy0, r0.quantile(0.25)['crown_sat'], r0.quantile(0.75)['crown_sat'], color='b', alpha=0.3)
# ax[1].plot(doy0, r0.median()['crown_sat'], 'b-', label='1990-2020, sparse')

# ax[1].fill_between(doy1, r1.quantile(0.25)['crown_sat'], r1.quantile(0.75)['crown_sat'], color='r', alpha=0.3)
# ax[1].plot(doy1, r1.median()['crown_sat'], 'r-', label='1990-2020, dense')

# ax[1].legend()
# ax[1].set_ylabel('crown saturation'); ax[1].set_xlabel('doy')


# ax[2].fill_between(doy0, r0.quantile(0.25)['litter_sat'], r0.quantile(0.75)['litter_sat'], color='b', alpha=0.3)
# ax[2].plot(doy0, r0.median()['litter_sat'], 'b-', label='1990-2020, sparse')

# ax[2].fill_between(doy1, r2.quantile(0.25)['litter_sat'], r2.quantile(0.75)['litter_sat'], color='r', alpha=0.3)
# ax[2].plot(doy1, r2.median()['litter_sat'], 'r-', label='2069-2099 RCP4.5, sparse')
# ax[2].legend()

# #%%

# fig, ax = plt.subplots(3,1)
# fig.set_size_inches(11.5, 11.5)


# ax[0].fill_between(doy0, f0.quantile(0.25)['Prec']*86400, f0.quantile(0.75)['Prec']*86400, color='b', alpha=0.3)
# ax[0].plot(doy0, f0.mean()['Prec']*86400, 'b-', label='1990-2020')
# ax[0].set_ylabel('Prec (mm d-1)')

# ax[0].fill_between(doy0, f1.quantile(0.25)['Prec']*86400, f1.quantile(0.75)['Prec']*86400, color='r', alpha=0.3)
# ax[0].plot(doy0, f1.mean()['Prec']*86400, 'r-', label='2069-2099 RCP4.5')
# ax[0].set_ylabel('Prec (mm d-1)')

# ax[0].legend()
# #plt.savefig(r'Data/Pallas_pine_climate.png', dpi=900))

# #%%
# res0['year'] = res0.index.year
# res0['month'] = res0.index.month
# res0['day'] = res0.index.day

# res1['year'] = res1.index.year
# res1['month'] = res1.index.month
# res1['day'] = res1.index.day

# r0 = res0.groupby('month')
# r1 = res1.groupby('month')
# #f0 = forc_present.groupby('month')
# #f1 = forc_future.groupby('month')

# doy0 = r0.mean().index.values
# doy1 = r1.mean().index.values

# fig, ax = plt.subplots(3,1)
# fig.set_size_inches(11.5, 11.5)

# ax[0].fill_between(doy0, r0.quantile(0.25)['litter_sat'], r0.quantile(0.75)['litter_sat'], color='b', alpha=0.3)
# ax[0].plot(doy0, r0.median()['litter_sat'], 'b-', label='1990-2020, sparse')

# ax[0].fill_between(doy1, r1.quantile(0.25)['litter_sat'], r1.quantile(0.75)['litter_sat'], color='r', alpha=0.3)
# ax[0].plot(doy1, r1.median()['litter_sat'], 'r-', label='1990-2020, dense')

# ax[0].legend()
# ax[0].set_ylabel('litter saturation'); ax[0].set_xlabel('doy')

# ax[1].fill_between(doy0, r0.quantile(0.25)['crown_sat'], r0.quantile(0.75)['crown_sat'], color='b', alpha=0.3)
# ax[1].plot(doy0, r0.median()['crown_sat'], 'b-', label='1990-2020, sparse')

# ax[1].fill_between(doy1, r1.quantile(0.25)['crown_sat'], r1.quantile(0.75)['crown_sat'], color='r', alpha=0.3)
# ax[1].plot(doy1, r1.median()['crown_sat'], 'r-', label='1990-2020, dense')

# ax[1].legend()
# ax[1].set_ylabel('crown saturation'); ax[1].set_xlabel('doy')


# ax[2].fill_between(doy0, r0.quantile(0.25)['swe'], r0.quantile(0.75)['swe'], color='b', alpha=0.3)
# ax[2].plot(doy0, r0.median()['swe'], 'b-', label='1990-2020, sparse')
# ax[2].fill_between(doy1, r1.quantile(0.25)['swe'], r1.quantile(0.75)['swe'], color='r', alpha=0.3)
# ax[2].plot(doy0, r1.median()['swe'], 'r-', label='1990-2020, dense')


