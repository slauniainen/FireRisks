# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 15:33:47 2022

module for reading and manipulating mNFI forest attributes (Mäkisara et al., 2016)

@author: Samuli Launiainen
"""
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import rasterio

EPS = np.finfo(float).eps

nodata = 32765.
#folder = r'D:\mVMI\2019'

#tiffs = glob.glob( r'D:/mVMI/2019/*.tif')

def sample_mNFI(maintype=1, sitetype=4, species=('p', 0.9)):
    
    files = {
            'maintype': r'D:/mVMI/2019\\paatyyppi_vmi1x_1519_M4.tif',
            'sitetype': r'D:/mVMI/2019\\kasvupaikka_vmi1x_1519_M4.tif',
            'age': r'D:/mVMI/2019\\ika_vmi1x_1519_M4.tif',
            'ba': r'D:/mVMI/2019\\ppa_vmi1x_1519_M4.tif',
            'd50': r'D:/mVMI/2019\\keskilapimitta_vmi1x_1519_M4.tif',
            'height': r'D:/mVMI/2019\\keskipituus_vmi1x_1519_M4.tif',
            'vol': r'D:/mVMI/2019\\tilavuus_vmi1x_1519_M4.tif',
            'vol_p': r'D:/mVMI/2019\\manty_vmi1x_1519_M4.tif',
            'vol_s': r'D:/mVMI/2019\\kuusi_vmi1x_1519_M4.tif',
            'fol_p': r'D:/mVMI/2019\\bm_manty_neulaset_vmi1x_1519_M4.tif',
            'fol_s': r'D:/mVMI/2019\\bm_kuusi_neulaset_vmi1x_1519_M4.tif',
            'fol_d': r'D:/mVMI/2019\\bm_lehtip_neulaset_vmi1x_1519_M4.tif'
            }
    
    raw = {k: [] for k in files.keys()}
    
    for k, fname in files.items():
        src = rasterio.open(fname)
        data = src.read(1).astype(float)
        src.close()
        raw[k] = data
    
    
    # resample
    if species[0] == 'pine':   
        f = raw['vol_p'] / (raw['vol'] + 1.0)
        raw['fol'] = raw['fol_p']
        raw['f_vol'] = f
        del raw['vol_s'], raw['fol_s'], raw['fol_d'], raw['fol_p']
        
    elif species[0] == 'spruce':
        f = raw['vol_s'] / (raw['vol'] + 1.0)
        raw['f_vol'] = f
        raw['fol'] = raw['fol_s']
        del raw['vol_p'],  raw['fol_s'], raw['fol_d'], raw['fol_p']
    else:
        f = (raw['vol'] - raw['vol_p'] - raw['vol_s']) /  (raw['vol'] + 1.0)
        raw['f_vol'] = f
        raw['fol'] = raw['fol_d']
        del raw['vol_p'], raw['vol_s'],  raw['fol_s'], raw['fol_d'], raw['fol_p']
        
    mask = np.where((raw['maintype'] == maintype) & (raw['sitetype'] == sitetype) & (f > species[1]))

    for k, v in raw.items():
        raw[k] = v[mask]
    
    raw['height'] = 0.1 * raw['height']
    
    return raw

def percentile_dbh_distr(species, d50, G, age, soil, do):
    """ 
    Kangas A. & Maltamo M. 2000. Silva Fennica
    percentile dbh distibutions for scots pine & norway spruce
    Args:
        species - 'pine', 'spruce', 'birch (ADD!)'
        d50 - stand median dbh (cm)
        G - stand basal area (m2 ha-1)
        age - stand age (years)
        Soil - dummy: 1 for mineral soils with fertility class >= mesic
    Returns:
        dbh -frequency (array, ha-1)
        ba - distribution (array, m2 ha-1 per size class)
    """
    
    def fun(species, d50, G, age, soil):
        """
        diameters at 0...100 percentiles of the dbh-distribution
        """
        
        params = {'pine':
              #interc., ln(d50),ln(age), ln(age/G), Soil*Ln(age), Soil, quantile
            [[-1.3776, 1.3656, -0.1886, 0.0300, 0., 0., 0.], #d0
             [-0.8903, 1.3328, -0.1405, 0., 0., -0.0281, 0.1], 
             [-0.5176, 1.1976, -0.0830, 0., 0, -0.0351, 0.2],
             [-0.3422, 1.1430, -0.0590, 0., 0., -0.0151, 0.3],
             [-0.1550, 1.0561, -0.0218, 0., 0., -0.0054, 0.4],
             [0., 1.0, 0., 0., 0., 0., 0.5],                 #d50 from inputs     
             [0.1576, 0.9264, 0.0297, 0., 0., 0., 0.6],
             [0.3282, 0.8452, 0.0621, 0., 0., 0., 0.7],
             [0.4505, 0.7881, 0.0896, 0., 0., 0., 0.8],
             [0.6887, 0.7214, 0.1025, 0., 0., 0., 0.9],
             [0.7921, 0.7032, 0.1045, 0., 0., 0., 0.95],
             [1.0369, 0.5889, 0.1668, -0.0215, -0.0042, 0., 1.0]], #d100
           'spruce':
              #interc., ln(d50),ln(age), ln(age/G), ln(G), Soil, quantile
            [[-0.3561, 0.8351, 0., -0.1178, -0.1261, 0.0, 0.0],
             [-0.2120, 0.8830, 0., -0.0736, 0., 0., 0.1],
             [-0.1667, 0.9679, 0., -0.0789, 0., 0., 0.2],
             [-0.3199, 1.0528, 0., -0.0379, 0., 0., 0.3],
             [-0.1315, 1.0266, 0., -0.0313, 0., 0., 0.4],
             [0., 1.0, 0., 0., 0., 0., 0],
             [0.1766, 0.9688, 0., 0., 0., 0., 0.6],
             [0.3237, 0.8964, 0.0348, 0., 0.,0., 0.7],
             [0.4768, 0.8381, 0.0603, 0., 0., 0., 0.8],
             [0.7771, 0.7502, 0.0792, 0., 0., -0.0392, 0.9],
             [0.9005, 0.7016, 0.1014, 0., 0., -0.0409, 0.95],
             [1.3823, 0.6241, 0.0832, 0., 0., -0.0682, 1.0]],
              
            }
            
        p = np.array(params[species])      
        r = len(p)

        q  = p[:,-1] # quantiles
        
        ln_d = np.ones(r) * np.NaN
        
        if species == 'pine':
            for k in range(r):
                ln_d[k] = p[k,0] +  np.log(d50) * p[k,1] + np.log(age) * p[k,2] \
                    + np.log(age/(G+EPS)) * p[k,3] + soil * np.log(age) * p[k,4] + soil * p[k,5]
        
        elif species == 'spruce':
            for k in range(r):
                ln_d[k] = p[k,0] +  np.log(d50) * p[k,1] + np.log(age) * p[k,2] \
                    + np.log(age/(G+EPS)) * p[k,3] + np.log(G) * p[k,4] + soil * p[k,5]
                    
        d = np.exp(ln_d)
        
        F = interp1d(d, q, kind='cubic', bounds_error=False, fill_value=(0.0, 1.0))
        
        return F
    
    # mean dbh in each class
    dbh = 0.5 * (do[0:-1] + do[1:])
    tree_ba = np.pi * (0.005 * dbh)**2

    N = np.zeros((len(dbh), len(d50)))
    Ntot = np.zeros(len(d50))
    
    
    for k in range(len(d50)):
        if d50[k] > 1.0:
            # get cumulative dbh-basal area function
            F = fun(species, d50[k], G[k], age[k], soil)
            # evaluate at do-intervals
            y = np.maximum(0., F(do[1:]) - F(do[0:-1])) * G[k]
            # trees ha-1 per size class
            ntrees = y  / tree_ba
            N[:,k] = ntrees
            Ntot[k] = sum(ntrees)
            
    return dbh, N, Ntot