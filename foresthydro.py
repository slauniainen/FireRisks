# -*- coding: utf-8 -*-
"""
Simple water balance model to estimate forest (canopy and litter layer) rainfall 
interception and wetness dynamics

Two cases illustrated in respective Notebooks:
    1) predicting litter and crown moisture pdf's for snow-free season
    2) approximating rainfall interception fraction per event size in forests of
    different densities

Samuli Launiainen 5.5.2022; uses same waterbudget model than ongoing moss Nfix
model.

"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from waterbalance import Canopywaterbudget

#: machine epsilon
EPS = np.finfo(float).eps

class Model(object):
    
    def __init__(self, sitepara, wpara):
        """
        Creates model object

        Returns:
            self
        """
        """
        water_p = {'interception': {'LAI': 4.0, 'wmax': 0.2, 'alpha': 1.28},
          'snowpack': {'Tm': 0.0, 'Km': 2.9e-5, 'Kf': 5.8e-6, 'R': 0.05, 
                       'Tmax': 0.5, 'Tmin': -0.5, 'Sliq': 0.0, 'Sice': 0.0},
          'organic_layer': {'DM': 0.1, 'Wmax': 10.0, 'Wmin': 0.1, 'Wcrit': 4.0,
                            'alpha': 1.28, 'Wliq': 10.0}                  
          }
         """
         
        self.dt = 86400.0 # s
        
        print(wpara, sitepara)
        wpara['interception']['LAI'] = sitepara['LAI']
        wpara['organic_layer']['DM'] = sitepara['DM_litter'] # kg m-2

        
        self.LAI = sitepara['LAI'] # (m2m-2)
        
        # waterbalance instance
        self.waterbalance = Canopywaterbudget(wpara)
        
    def run(self, forc):
        """
        Solves water budget at daily timestep
        Args:
            forc - weather forcing
        Returns:
            throughfall
            canopy water content
            litter water content
            snow water equivalent
        """
        
        # solve water balance
        
        # fluxes kg m-2(ground) s-1
        Trfall, Ecan, Ebl = self.waterbalance.run(self.dt, forc)
        canopy_water_storage = self.waterbalance.interception.W
        litter_water_content = self.waterbalance.bottomlayer.Wliq
        swe = self.waterbalance.snow.SWE
    
        return Trfall, canopy_water_storage, litter_water_content, swe
    
    def set_state(self, LAI=None, litter_dm=None):
        """
        This changes forest LAI [m2 m-2] and organic layer dry mass [kg m-2]
        """
        if LAI:
            self.waterbalance.manage_forest(LAI=LAI)
        if litter_dm:
            self.waterbalance.manage_forest(DM=litter_dm)
            

def driver(forc, sitepara, water_p):
    
    task = Model(sitepara, water_p)
    
    N = len(forc)
    
    results = {'date': forc['date'], 'Prec': np.ones((N,1))*np.NaN, 'Trfall': np.ones((N,1))*np.NaN,
               'canopy_water_storage': np.ones((N,1))*np.NaN, 
               'litter_water_content': np.ones((N,1))*np.NaN, 'swe': np.ones((N,1))*np.NaN}
    
    for k in range(N):
        print(k/N*100)
        f = forc.iloc[k]
        
        if f['doy'] == 1:
            task.set_state()
        
        trfall, cws, lmois, swe = task.run(forc.iloc[k].to_dict())
        
        results['Trfall'][k] = trfall
        results['canopy_water_storage'][k] = cws
        results['litter_water_content'][k] = lmois
        results['swe'][k] = swe
    
        results['Prec'][k] = forc['Prec'].iloc[k]
    
    
    for k in results.keys():
        results[k] = results[k].ravel()
    results = pd.DataFrame.from_dict(results)
    results.index = forc.index
    #results = results.drop(columns={'date'})
    results['doy'] = forc.doy
    return results, task

        
        
        
        
        