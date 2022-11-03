# -*- coding: utf-8 -*-
"""
Created on Thu May  5 15:45:51 2022

@author: Samuli Launiainen

Wind speed profile and aerodynamic + boundary-layer conductances in canopies.
"""
import numpy as np
       
def aerodynamics(LAI, hc, Uo, w=0.01, zm=2.0, zg=0.5, zos=0.01):
    """
    computes wind speed at ground and canopy + boundary layer conductances
    Computes wind speed at ground height assuming logarithmic profile above and
    exponential within canopy
    Args:
        LAI - one-sided leaf-area /plant area index (m2m-2)
        hc - canopy height (m)
        Uo - mean wind speed at height zm (ms-1)
        w - leaf length scale (m)
        zm - wind speed measurement height above canopy (m)
        zg - height above ground where Ug is computed (m)
        zos - forest floor roughness length, ~ 0.1*roughness element height (m)
    Returns:
        ra - canopy aerodynamic resistance (sm-1)
        rb - canopy boundary layer resistance (sm-1)
        ras - forest floor aerod. resistance (sm-1)
        ustar - friction velocity (ms-1)
        Uh - wind speed at hc (ms-1)
        Ug - wind speed at zg (ms-1)
    SOURCE:
       Cammalleri et al. 2010 Hydrol. Earth Syst. Sci
       Massman 1987, BLM 40, 179 - 197.
       Magnani et al. 1998 Plant Cell Env.
    """
    zm = hc + zm  # m
    zg = np.minimum(zg, 0.1 * hc)
    kv = 0.4  # von Karman constant (-)
    beta = 285.0  # s/m, from Campbell & Norman eq. (7.33) x 42.0 molm-3
    alpha = LAI / 2.0  # wind attenuation coeff (Yi, 2008 eq. 23)
    d = 0.66*hc  # m
    zom = 0.123*hc  # m
    zov = 0.1*zom
    zosv = 0.1*zos

    # solve ustar and U(hc) from log-profile above canopy
    ustar = Uo * kv / np.log((zm - d) / zom) 
    Uh = ustar / kv * np.log((hc - d) / zom)
    
    # U(zg) from exponential wind profile
    zn = np.minimum(zg / hc, 1.0)  # zground can't be above canopy top
    Ug = Uh * np.exp(alpha*(zn - 1.0))

    # canopy aerodynamic & boundary-layer resistances (sm-1). Magnani et al. 1998 PCE eq. B1 & B5
    #ra = 1. / (kv*ustar) * np.log((zm - d) / zom)
    ra = 1./(kv**2.0 * Uo) * np.log((zm - d) / zom) * np.log((zm - d) / zov)    
    rb = 1. / LAI * beta * ((w / Uh)*(alpha / (1.0 - np.exp(-alpha / 2.0))))**0.5

    # soil aerodynamic resistance (sm-1)
    ras = 1. / (kv**2.0*Ug) * (np.log(zg / zos))*np.log(zg / (zosv))
    
    #print('ra', ra, 'rb', rb)
    ra = ra + rb
    return ra, rb, ras, ustar, Uh, Ug


def wind_profile(LAI, hc, Uo, z, zm=2.0):
    """
    Computes wind speed at ground height assuming logarithmic profile above and
    hyperbolic cosine profile within canopy
    INPUT:
        LAI - one-sided leaf-area /plant area index (m2m-2)
        hc - canopy height (m)
        Uo - mean wind speed at height zm (ms-1)
        zm - wind speed measurement height above canopy (m)
    OUTPUT:
        U - wind speed profile(ms-1)
        ustar - friction velocity (ms-1)
        Uh - wind speed at canopy top (ms-1)
    SOURCE:
       Cammalleri et al. 2010 Hydrol. Earth Syst. Sci
       Massman 1987, BLM 40, 179 - 197.
    """

    k = 0.4  # von Karman const
    Cd = 0.2  # drag coeff
    alpha = 1.5  # (-)

    zm = zm + hc
    d = 0.66*hc
    zom = 0.123*hc
    beta = 4.0 * Cd * LAI / (k**2.0*alpha**2.0)
    # solve ustar and U(hc) from log-profile above canopy
    ustar = Uo * k / np.log((zm - d) / zom)  # m/s

    U = np.ones(len(z))*np.NaN

    # above canopy top wind profile is logarithmic
    U[z >= hc] = ustar / k * np.log((z[z >= hc] - d) / zom)

    # at canopy top, match log and exponential profiles
    Uh = ustar / k * np.log((hc - d) / zom)  # m/s

    # within canopy hyperbolic cosine profile
    U[z <= hc] = Uh * (np.cosh(beta * z[z <= hc] / hc) / np.cosh(beta))**0.5

    return U, ustar, Uh 