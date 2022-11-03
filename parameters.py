# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 13:53:57 2021

@author: 03081268
"""
import pandas as pd
import os
import matplotlib.pyplot as plt

def read_sitedata(ffile, siteid=None):
    """
    returns site data as dataframe
    """
    
    dat = pd.read_csv(ffile, sep=';', header='infer')
    dat.index = dat['ID']
    
    if siteid:
        dat = dat.iloc[siteid]
        
    return dat

def read_forcing(scenario, grid_id, fdate=None, ldate=None):
    """
    returns meteorological data as dataframe
    """
    
    ffile =  os.path.join(r'Data\forcing', 'CanESM2_' + scenario, 'weather_id_' +str(grid_id) + '.csv')
    print(ffile)
    dat = pd.read_csv(ffile, sep=',', header='infer')
    dat = dat.rename(columns={'TAir': 'Tair', 'Precip': 'Prec', 'global_radiation': 'Rg'})
    dat.Prec  /= 86400.

    if fdate:
        dat = dat[dat['date']>=fdate]
    if ldate:
        dat = dat[dat['date']<=ldate]
        
    return dat


#%% parameters
sfile = r'Data/ICP_L2_sitedata.csv'
rcp_scen = 'rcp45'
#grid_id = 3812 #[68.6N / 25.2E]
#grid_id = 3787 #[65.6N / 25.6E]
#grid_id = 3286 #[63.5N / 24.4E]
grid_id = 3269 #[61.8N / 24.4E]


sitedata = read_sitedata(sfile, siteid=3)
forc = read_forcing(rcp_scen, grid_id, fdate='1990-01-01', ldate='2020-12-31')

# canopy waterbudget
water_p = {'interception': {'wmax': 0.2, 'alpha': 1.28},
          'snowpack': {'Tm': 0.0, 'Km': 2.9e-5, 'Kf': 5.8e-6, 'R': 0.05, 
                       'Tmax': 0.5, 'Tmin': -0.5, 'Sliq': 0.0, 'Sice': 0.0},
          'organic_layer': {'Wmax': 10.0, 'Wmin': 0.1, 'Wcrit': 4.0,
                            'alpha': 1.28, 'Wliq': 10.0}                  
          }

nitrogen_p = { 'Hylo': {'potrate': [25.6, 12.8], 'DM': [sitedata['Hylo_u'], sitedata['Hylo_l']], 'wresponse': [0.44, 7.0]},
               'Pleu': {'potrate':[27.0, 7.3], 'DM': [sitedata['Pleu_u'], sitedata['Pleu_l']], 'wresponse': [0.44, 7.0]},
               'Dicr': {'potrate': [0.2, 0.0], 'DM': [sitedata['Dicr_u'], sitedata['Dicr_l']], 'wresponse': [0.44, 7.0]},
             }


# bucket
#soil_p = {'depth': 0.5, 'Ksat': 1e-5,
#          'pF': {'ThetaS': 0.50, 'ThetaR': 0.03, 'alpha': 0.06, 'n': 1.35},# Hyytiala A-horizon pF curve 
#                 'MaxPond': 0.0, 'Wliq': 0.4
#         }