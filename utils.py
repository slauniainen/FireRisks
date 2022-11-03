# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 13:53:57 2021

@author: 03081268
"""
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import sys

EPS = np.finfo(float).eps

def read_sitedata(ffile, siteid=None):
    """
    returns site data as dataframe. Biomasses are as kg ha-1 dm
    """
    
    dat = pd.read_csv(ffile, sep=';', header='infer')
    #dat.index = dat['ID']
    
    if siteid:
        dat = dat[dat.ID==siteid]
        dat = dat.squeeze()
        
    return dat

def read_climate_forcing(scenario, grid_id, fdate=None, ldate=None):
    """
    returns meteorological data as dataframe
    """
    
    ffile =  os.path.join(r'Data\forcing', 'CanESM2_' + scenario, 'weather_id_' +str(grid_id) + '.csv')
    print(ffile)
    dat = pd.read_csv(ffile, sep=',', header='infer')
    dat = dat.rename(columns={'TAir': 'T', 'Precip': 'Prec', 'global_radiation': 'Rg'})
    dat.Prec  /= 86400.

    if fdate:
        dat = dat[dat['date']>=fdate]
    if ldate:
        dat = dat[dat['date']<=ldate]

    tvec = pd.to_datetime(dat['date'], yearfirst=True)
    dat.index = tvec
    dat.date = tvec
    return dat

def read_forcing(forc_fp, start_time, end_time,
                 dt=1800.0, na_values='NaN', sep=';'):
    """
    Reads forcing data from to dataframe
    Args:
        forc_fp (str): forcing file name
        start_time (str): starting time [yyyy-mm-dd], if None first date in
            file used
        end_time (str): ending time [yyyy-mm-dd], if None last date
            in file used
        dt (float): time step [s], if given checks
            that dt in file is equal to this
        na_values (str/float): nan value representation in file
        sep (str): field separator
    Returns:
        Forc (dataframe): dataframe with datetime as index and cols read from file
    """


    dat = pd.read_csv(forc_fp, header='infer', na_values=na_values, sep=sep)

    # set to dataframe index
    tvec = pd.to_datetime(dat[['year', 'month', 'day', 'hour', 'minute']])
    tvec = pd.DatetimeIndex(tvec)
    dat.index = tvec

    dat = dat[(dat.index >= start_time) & (dat.index <= end_time)]

    # convert: H2O mmol / mol --> mol / mol; Prec kg m-2 in dt --> kg m-2 s-1
    dat['H2O'] = 1e-3 * dat['H2O']
    dat['Prec'] = dat['Prec'] / dt
    
    es, _ = e_sat(dat['Tair'])
    dat['VPD'] = 1e-3 * (es - dat['H2O'] * dat['P']) # kPa
    dat['VPD'] = np.maximum(EPS, dat['VPD'])
    dat['RH'] = 100 * dat['H2O'] * dat['P'] / es
    dat['RH'] = np.maximum(EPS, dat['RH'])
                    
    cols = ['doy', 'Prec', 'P', 'Tair', 'Tdaily', 'U', 'Ustar', 'H2O', 'VPD', 'RH', 'CO2', 'Zen',
            'LWin', 'diffPar', 'dirPar', 'diffNir', 'dirNir']
    # these for phenology model initialization
    if 'X' in dat:
        cols.append('X')
    if 'DDsum' in dat:
        cols.append('DDsum')

    # for case for bypassing soil computations
    if 'Tsh' in dat:
        cols.append('Tsh')
    if 'Wh' in dat:
        cols.append('Wh')
    if 'Tsa' in dat:
        cols.append('Tsa')
    if 'Wa' in dat:
        cols.append('Wa')
    if 'Rew' in dat:
        cols.append('Rew')        
    # Forc dataframe from specified columns
    Forc = dat[cols].copy()

    # Check time step if specified
    if len(set(Forc.index[1:]-Forc.index[:-1])) > 1:
        sys.exit("Forcing file does not have constant time step")
    if (Forc.index[1] - Forc.index[0]).total_seconds() != dt:
        sys.exit("Forcing file time step differs from dt given in general parameters")

    return Forc

def read_data(ffile, start_time=None, end_time=None, na_values='NaN', sep=';'):

    dat = pd.read_csv(ffile, header='infer', na_values=na_values, sep=sep)
    # set to dataframe index
    tvec = pd.to_datetime(dat[['year', 'month', 'day', 'hour', 'minute']])
    tvec = pd.DatetimeIndex(tvec)
    dat.index = tvec

    # select time period
    if start_time == None:
        start_time = dat.index[0]
    if end_time == None:
        end_time = dat.index[-1]

    dat = dat[(dat.index >= start_time) & (dat.index <= end_time)]

    return dat

def e_sat(T):
    """
    Computes saturation vapor pressure (Pa), slope of vapor pressure curve
    [Pa K-1]  and psychrometric constant [Pa K-1]
    IN:
        T - air temperature (degC)
    OUT:
        esa - saturation vapor pressure in Pa
        s - slope of saturation vapor pressure curve (Pa K-1)
    SOURCE:
        Campbell & Norman, 1998. Introduction to Environmental Biophysics. (p.41)
    """

    esa = 611.0 * np.exp((17.502 * T) / (T + 240.97))  # Pa
    s = 17.502 * 240.97 * esa / ((240.97 + T)**2)

    return esa, s