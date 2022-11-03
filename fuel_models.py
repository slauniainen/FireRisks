# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 15:11:14 2022

Module of dead fuel moisture models

@author: Samuli Launiainen
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

EPS = np.finfo(float).eps


def FM_deDios(D, fm0, fm1, tau):
    """
    deDios et al. 2015 Agric. For. Met.
    Args:
        D - daily maximum VPD (kPa), scalar or array
        fm0 - minimum fuel moisture (% or g g-1)
        fm0 - maximum fuel moisture (% or g g-1)
        tau - time constant (kPa-1)
    Returns:
        fm (% or g g-1, eq. 1 or 7)
    """
    
    fm = fm0 + fm1 * np.exp(-tau * D)
    return fm


def FM_equilibration(T, RH, Prec, dt, tau, fm_max=None, Prec_lim=1.0):
    """
    Based on Catchpole et al. 2001 Int. J. Wildl. Fire &
        Matthews 2014 Int. J. Wildl. Fire with heuristic impact of precipitation
    Args:
        T - air temperature (degC), array
        RH - air relative humidity (%), array
        Prec - precipitation rate (mm/h)
        dt - timestep (hours)
        tau - time constant (hours)
        Prec_lim - threshold precipitation rate (mm/h) to saturate fuel 
    Returns:
        fm (units)
        EMC - equilibrium moisture content (units)
    """
    
    #R = 8.314 # J mol-1 K-1
    #Mw = 18.015e-3 # molar mass of water (kg mol-1)
    
    # the papers seem to give parameters in units that need R and Mw to be
    R = 1.987 # cal mol-1
    Mw = 18.015 # g mol-1
    #... do the unit conversion later on
    
    # parameters a & b (g g-1) are based on adsoption/desorption isotherms; 
    # fuel specific. See Nelson 2001; Blackmarr 1971; Anderson 1990 https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=1055&context=govdocs_forest
    # do we have data for boreal forest litter/moss layer?
   
    # these are from my generic fit
    a = 0.32 # g/g range 0.26...0.33, Catchpole et al. 2001
    b = -0.062 # (g/g), range -0.05...0.08, Catchpole et al. 2001, Nelson 2000    
    
    # equilibrium moisture content (g g-1); Catchpole eq. 9.
    # In reality, a and b should be fuel-specific --> adjust to determine range of FM
    
    EMC = a + b*np.log(- R*(T + 273.15) / Mw * np.log(RH/100.0))

    N = len(T)
    
    fm = np.zeros(N)
    fm[0] = EMC[0] # assume initial state equal to EMC
    
    # Stetson-Harrison: maximum water content would correspond to EMC at 25degC and 99%RH
    if not fm_max:
        fm_max = a + b*np.log(- R*(25.0 + 273.15) / Mw * np.log(99.0 / 100.0))
    
    # compute fm dynamics
    for k in range(1, N):
        # heuristic effect of precipitation
        if Prec[k] > Prec_lim:
            #fm[k] = fm_max
            fm[k] = np.max([fm[k-1], fm_max * np.min([Prec[k] / Prec_lim, 1.0])])
        
        else:
            # Catchpole eq. 7; in principle same as below
            # L = np.exp(- dt / (2 * tau))
            #fm[k] = L**2 * fm[k-1] + L*(1-L)*EMC[k-1] + (1-L)*EMC[k]
        
            fm[k] = EMC[k] + (fm[k-1] - EMC[k])*np.exp(-dt/tau)
        
    #-- for testing
    fig, ax = plt.subplots(2,1)
    ax[0].set_title('euilibration')
    ax[0].plot(T, 'r-'); ax[1].set_ylabel('T (degC)')
    ax0b = ax[0].twinx()
    ax0b.plot(RH, 'k-'); ax0b.set_ylabel('RH (%)')
    ax[1].plot(fm, 'r-', label='fm, tau=%.1f' %tau)
    ax[1].plot(EMC, 'k-', label='EMC')
    ax[1].set_ylabel('fm (g g-1)')
    #ax[2].plot(c)
    return fm, EMC

def FM_catchpole_modified(T, RH, Prec, dt, tau, Prec_lim):
    """
    Catchpole et al. 2001 Int. J. Wildl. Fire
    Testing to include rainfall-effects (instantaneous wetting to saturation)
    Args:
        T - air temperature (degC), array
        RH - air relative humidity (%), array
        Prec - precipitation rate (mm h-1), array
        dt - timestep (hours)
        tau - time constant (hours)
        Prec_lim - threshold precipitation to saturate the fuel
    Returns:
        fm (units)
        EMC - equilibrium moisture content (units)
    """
    
    #R = 8.314 # J mol-1 K-1
    #Mw = 18.015e-3 # molar mass of water (kg mol-1)
    
    
    # the papers seem to give parameters in units that need R and Mw to be
    R = 1.987 # cal mol-1
    Mw = 18.015 # g mol-1
    #... do the unit conversion later on
    
    # parameters a & b (g g-1) are based on adsoption/desorption isotherms; 
    # fuel specific. See Nelson 2001; Blackmarr 1971; Anderson 1990 https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=1055&context=govdocs_forest
    
    # do we have data for boreal forest litter/moss layer?
    a = 0.30 # g/g range 0.26...0.33, Catchpole et al. 2001
    b = -0.06 # (g/g), range -0.05...0.08, Catchpole et al. 2001, Nelson 2000    
    
    # handle inputs; restrict RH to range [1, 99 %].
    #RH = np.minimum(np.maximum(1.0, RH), 99.0)
    
    # equilibrium moisture content (g g-1); Catchpole eq. 9. a and b should be fuel-specific
    EMC = a + b*np.log(- R*(T + 273.15) / Mw * np.log(RH/100.0))
    
    fm_max = np.max(a + b*np.log(- R*(T + 273.15) / Mw * np.log(99.0 / 100.0)))
    
    N = len(T)
    
    fm = np.zeros(N)
    fm[0] = EMC[0] # assume initial state equal to EMC
    
    L = np.exp(- dt / (2 * tau))

    for k in range(1, N):
        
        #...handwaving: if Prec > threshold, set fm to max; else proportion to max.
        if Prec[k] > Prec_lim:
            fm[k] = fm_max
            #fm[k] = np.max([fm[k-1], fm_max * np.min([Prec[k] / Prec_lim, 1.0])])
        else:
            # Catchpole eq. 7
            fm[k] = L**2 * fm[k-1] + L*(1-L)*EMC[k-1] + (1-L)*EMC[k]

    #-- for testing
    
    fig, ax = plt.subplots(2,1)
    ax[0].set_title('catchpole + rain')
    ax[0].plot(T, 'r-'); ax[1].set_ylabel('T (degC)')
    ax0b = ax[0].twinx()
    ax0b.plot(RH, 'k-'); ax0b.set_ylabel('RH (%)')
    ax[1].plot(fm, 'r-', label='fm, tau=%.1f' %tau)
    ax[1].plot(EMC, 'k-', label='EMC')
    ax[1].set_ylabel('fm (g g-1)')
    
    return fm, EMC
