# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 09:54:34 2021

@author: 03081268
"""
import matplotlib.pyplot as plt
import numpy as np

#: machine epsilon
EPS = np.finfo(float).eps
#: [kg m\ :sup:`-3`\ ], water density
WATER_DENSITY = 1.0e3

class Canopywaterbudget(object):
    """
    Combines Interception, Snowpack, Organiclayer, CanopyConductance and Rootzone
    """
    
    def __init__(self, dt, LAI, params):
        
        self.dt = dt

        self.interception = Interception(LAI, params['interception'])
        self.snow = Snowpack(params['snowpack'])
        self.bottomlayer = OrganicLayer(params['organic_layer'])
    
        self.canopy = CanopyTranspiration(LAI, params['transpiration'])
        self.soil = Bucket(params['soil'])
        
    def run(self, forcing):
        """
        solves sequantial canopy and soil bucket water balances
        Args:
            dt - timestep (s)
            doy - day of year
            T - air temperature (degC)
            D - vapor pressure deficit (kPa)
            Prec - precipitation rate (mm s-1)
            Rg - global radiation (W m-2)
        Returns:
            Prec - throughfall to soil (mm s-1)
            Ecan - canopy evaporation rate (mm s-1)
            Ebl - bottom layer evaporation rate (mm s-1)
            Tr - transpiration rate (mm s-1)
            litter moisture - bottom layer water content (g g-1)
            soil water content - root zone (m3 m-3)
            soil water potential - root zone (Mpa)
        """

        doy = forcing['doy']
        T = forcing['T']
        Prec = forcing['Prec']
        D = forcing['VPD']
        Rg = forcing['Rg']
        Qp = 0.5 * Rg # half of global radiation is PAR
        
        # approximate net radiation at canopy and ground layer
        f = np.exp(-0.6 * self.interception.LAI)
        Rnc = 0.8 * (1-f) * Rg # canopy layer
        Rng = 0.8 * f * Rg # ground layer
        
        # interception
        Prec, Interc, Ecan, imbe = self.interception.run(self.dt, T, Prec, Rnc)
        
        # snowpack
        Prec, mbe = self.snow.run(self.dt, T, Prec)
        
        # bottom layer
        Prec, Ebl = self.bottomlayer.run(self.dt, Prec, T, Rng, Snowcover=self.snow.SWE)
        
        # compute canopy transpiration rate (m s-1)
        Tr = self.canopy.transpiration_rate(doy, Qp, T, D, Rew=self.soil.Rew) 
        
        #... that cannot exceed available soil water
        Tr = np.minimum(Tr, WATER_DENSITY * self.soil.WatSto / self.dt)
        
        # solve root zone water balance
        Infil, Roff, Drain, mbe_bu = self.soil.run(self.dt, rr=Prec, et=Tr, latflow=0.0)
        
        #print(Prec, Infil, mbe_bu)
        
        # precipitation, canopy evaporation rate, litter layer evaporation rate, transpiration rate,
        # relative canopy storage, litter water content, soil water content, water potential
        
        return Prec, Ecan, Ebl, Tr, self.interception.W / self.interception.Wmax,\
                self.bottomlayer.Wliq, self.soil.Wliq, self.soil.h
    
    def manage_forest(self, LAI=None, DM=None):
        """
        adjust LAI and/or bottom layer dry mass (DM)
        Args:
            LAI - float
            DM - float
        """
        
        if LAI:
            self.canopy.LAI = LAI
            self.interception.LAI = LAI
            self.interception.Wmax = self.interception.wmax * LAI
            self.interception.W = np.minimum(self.interception.W, self.interception.Wmax)
        if DM:
            self.bottomlayer.DM = DM
            self.bottomlayer.WatSto = self.bottomlayer.Wliq * self.bottomlayer.DM

# submodels
class Interception():
    """
    Big-leaf model for rainfall interception in vegetation canopy
    Evaporation follows Priestley-Taylor equilibrium evaporation
    """
    def __init__(self, LAI, p):
        self.wmax =p['wmax'] # max storage (mm/LAI)
        self.LAI = LAI
        self.alpha = p['alpha'] # priestley-taylor alpha (-)
        
        self.Wmax = self.wmax * self.LAI
        self.W = 0.0  # storage mm
        
    def run(self, dt, T, Prec, AE):
        """
        Calculates canopy rainfall water interception and canopy water storage changes
        Args: 
            self - object
            T - air temperature (degC)
            Prec - precipitation rate (mm s-1 = kg m-2 s-1)
            AE - available energy (Wm-2)
        Returns:
            self - updated state W
            Trfall - thoughfall rate to forest floor (mm s-1)
            Interc - interception rate (mm s-1)
            Evap - evaporation from canopy store (mm s-1)
            MBE - mass balance error (mm)
        """
        Wo = self.W  # initial storage                
        Prec = Prec * dt  # mm

        # interception (mm)
        Interc = (self.Wmax - self.W) * (1.0 - np.exp(-(1.0 / self.Wmax) * Prec))
        # new canopy storage, mm
        self.W = self.W + Interc
        # Throughfall to field layer mm
        Trfall = Prec - Interc

        # evaporate from canopy storage mm s-1
        if Prec > 0:
            erate = 0.0  # zero during precipitation events
        else:
            erate, L = eq_evap(AE, T)
            erate = np.maximum(0.0, self.alpha * erate / L)
        
        Evap = np.minimum(erate * dt, self.W)  # mm
        
        self.W = self.W - Evap  # new storage

        # canopy storage mass-balance error
        MBE = (self.W - Wo) - (Prec - Evap - Trfall)

        return Trfall / dt, Interc / dt, Evap / dt,  MBE

class Snowpack(object):
    """
    degree-day snow model
    """
    
    def __init__(self, p):
        """
        Args:
            p (dict)
        """
        
        self.Tm = p['Tm'] # melt temp, degC
        self.Km = p['Km'] # melt coeff, mm d-1 degC-1
        self.Kf = p['Kf'] # freeze coeff
        self.R = p['R']  # max fraction of liquid water in snow
        self.Tmin = p['Tmin'] # below all Prec is snow
        self.Tmax = p['Tmax'] # above all Prec is liquid
        
        # water storages (liquid and ice) kgm-2 = mm
        self.liq = p['Sliq']
        self.ice = p['Sice']
        self.SWE = self.liq + self.ice
        
    def run(self, dt, T, Prec):
        """ 
        run snowpack
        Args:
            dt - [s]
            T - air temp [degC]
            Prec - precipitation rate [kg m-2 s-1 = mm s-1]
        """
        Prec = Prec * dt
        SWE0 = self.SWE
        
        # state of Prec depends on T
        fliq = np.maximum(0.0, 
                          np.minimum(1.0, (T - self.Tmin) / (self.Tmax - self.Tmin)))

        if T >= self.Tm:
            melt = np.minimum(self.ice, self.Km * dt * (T - self.Tm))  # mm
            freeze = 0.0
        else:
            melt = 0.0
            freeze = np.minimum(self.liq, self.Kf * dt * (self.Tm - T))  # mm
            
        # amount of water as ice and liquid in snowpack
        ice = np.maximum(0.0, self.ice + (1 - fliq) * Prec + freeze - melt)
        liq = np.maximum(0.0, self.liq + fliq * Prec - freeze + melt)

        trfall = np.maximum(0.0, liq - ice * self.R)  # mm

        # new state
        self.liq = np.maximum(0.0, liq - trfall)  # mm, liquid water in snow
        self.ice = ice
        self.SWE = self.liq + self.ice
        
        # mass-balance error mm
        mbe = (self.SWE - SWE0) - (Prec - trfall)
        
        return trfall / dt, mbe
        
class OrganicLayer(object):
    """
    Primitive model for organic layer moisture budget
    """
    def __init__(self, p):
        self.DM = p['DM'] # kg DMm-2
        self.Wmax = p['Wmax'] # gH2O/gDM, this equals 'field capacity'
        self.Wmin = p['Wmin']
        self.Wcrit = p['Wcrit'] # evaporation rate starts to decr. from here
        self.alpha = p['alpha'] # priestley-taylor alpha
        
        # initial state
        self.Wliq = self.Wmax * p['Wliq'] # g g-1
        self.WatSto = self.Wliq * self.DM # kg m-2 = mm
        
    def run(self, dt, Prec, T, Rnet, Snowcover=0.0):
        """
        Args:
            dt (s)
            Prec (kg m-2 s-1)
            T (degC)
            Rnet (Wm-2)
            Snowcover - snow water equivalent, >0 sets evap to zero
        Returns:
            evap (kg m-2 s-1)
            trall (kg m-2 s-1)
        """
        
        # rates
        if Snowcover > 0 or Prec > 0:
            evap = 0.0
        else:
            f = np.minimum(1.0, self.Wliq / self.Wcrit) # relative evaporation rate decreases in dry conditions
            Eo, L = eq_evap(Rnet, T) # Wm-2
            evap = f * self.alpha * Eo / L # kg m-2 s-1
            evap = np.minimum((self.Wliq - self.Wmin) * self.DM / dt, evap)
            
        interc = np.minimum(Prec, (self.Wmax - self.Wliq) * self.DM / dt)
        
        trfall = Prec - interc
        
        self.WatSto += (interc - evap) * dt
        self.Wliq = self.WatSto / self.DM
        
        return trfall, evap

class CanopyTranspiration(object):
    """
    Canopy transpiration rate
    Canopy conductance model based on Launiainen et al. 2019 HESS

    Assumes: LAI = const. during the year and canopy well-coupled to the atmosphere
    
    """
    def __init__(self, LAI, p):
        """
        Args:
            LAI (m2m-2) - one-sided leaf-area index
            p - parameters(dict)
                light and LAI-response: 
                    Amax (umol m-2 s-1), light-saturated leaf photosynthetic rate
                    g1 (kPa^0.5), USO model parameter
                    kp (-), attenuation coefficient for PAR
                    q50 (Wm-2), half-sat. of leaf light PAR response
                response to relative plant available water: 
                    rw - canopy conductance starts to drop when Rew <= rw
                    rwmin - minimum relative canopy conductance when Rew = 0.0
                seasonal cycle (delayed spring recovery) 
                    smax    
                    tau
                    xo
                    fmin - minimum relative canopy conductance (wintertime)
        """
        self.LAI = LAI
        self.p = p
        
        self.X = 0.0
        self.fPheno = 1.0
        self.DDsum = 0.0
        
    def photoacclim(self, T):
        """
        computes new stage of temperature acclimation and phenology modifier of
        stomatal conductance
        Peltoniemi et al. 2015 Bor.Env.Res.
        Args:
            T = daily mean air temperature
        
        """

        self.X = self.X + 1.0 / self.p['tau'] * (T - self.X)  # degC
        S = np.maximum(self.X - self.p['xo'], 0.0)
        fPheno = np.maximum(self.p['fmin'],
                            np.minimum(S / self.p['smax'], 1.0))
        
        self.fPheno = fPheno
        
        return fPheno
    
    def degreeDays(self, T, doy):
        """
        Calculates and updates degree-day sum from the current mean Tair.
        INPUT:
            T - daily mean temperature (degC)
            doy - day of year 1...366 (integer)
        """
        To = 5.0  # threshold temperature
        if doy == 1:  # reset in the beginning of the year
            self.DDsum = 0.
        else:
            self.DDsum += np.maximum(0.0, T - To)
    
    def Gc(self, Qp, T, D, fPheno=1.0, Rew=1.0, CO2=380.0):
        """
        Canopy conductance, follows Launiainen et al. 2019 HESS
        Args:
            p - parameters (dict)
            LAI - 1-sided leaf-area index (m2m-2)
            Qp - incoming PAR (W m-2)
            Ta - air temperature (degC)
            D - vapor pressure deficit (kPa)
            fPheno - seasonal cycle modifier (-)
            Rew - relatively extractable soil water (-)
            CO2 - ambient CO2 (ppm)
        Returns:
            Gc - canopy conductance (canopy-integrated stomatal conductance)  (m s-1)
        """
 
        D = np.maximum(D, EPS) # VPD alwas >= 0
        
        rhoa = 101300.0 / (8.31 * (T + 273.15)) # mol m-3
        
        Amax = self.p['Amax'] # light-saturated leaf photosynthetic rate (umol m-2 s-1)
        g1 = self.p['g1'] # USO model parameter (unit?)
        kp = self.p['kp']  # (-) attenuation coefficient for PAR
        q50 = self.p['q50']  # Wm-2, half-sat. of leaf light PAR response
        rw = self.p['rw']  # rew parameter; start of conductance decay
        rwmin = self.p['rwmin']  # rew parameter; minimum relative conductance
    
        LAI = self.LAI
        
        # leaf level light-saturated stomatal conductance gs (m/s)
        gs = 1.6*(1.0 + g1 / np.sqrt(D)) * Amax / CO2 / rhoa
    
        """--- environmental modifiers """
    
        # fQ (canopy-level light response): Saugier & Katerji, 1991 Agric. For. Met., eq. 4. Leaf light response = Qp / (Qp + q50)
        fQ = 1./ kp * np.log((kp*Qp + q50) / (kp*Qp*np.exp(-kp * LAI) + q50 + EPS) )
    
        # the next formulation is from Leuning et al., 2008 WRR for daily Gc; they refer to 
        # Kelliher et al. 1995 AFM but the resulting equation is not exact integral of K95.        
        # fQ = 1./ kp * np.log((Qp + q50) / (Qp*np.exp(-kp*self.LAI) + q50))
    
        # soil moisture response
        fRew = np.minimum(1.0, np.maximum(Rew / rw, rwmin))
    
        # CO2 -response of canopy conductance, derived from APES-simulations
        # (Launiainen et al. 2016, Global Change Biology). relative to 380 ppm
        fCO2 = 1.0 - 0.387 * np.log(CO2 / 380.0)
        
        # canopy conductance (ms-1)
        Gc = gs * fQ * fRew * fCO2 * fPheno
        
        return Gc
    
    def transpiration_rate(self, doy, Qp, T, D, Rew=1.0, Pamb=101.13):
        """
        canopy transpiration rate (m s-1)
        """
        # note -these are daily functions!
        self.degreeDays(T, doy)
        self.photoacclim(T)
        
        # canopy conductance
        Gc= self.Gc(Qp, T, D, Rew=Rew)
 
        Tr = Gc * D / Pamb
        
        return Tr
                
class Bucket(object):
    """
    Single-layer soil water bucket model (loosely following Guswa et al, 2002 WRR)
    """

    def __init__(self, p):
        """
        Args
            D - depth [m]
            Ksat - hydr. cond. [m/s]
            poros - porosity [vol/vol]
            pF - vanGenuchten water retention parameters {dict}
            MaxPond - maximum ponding depth above bucket [m]
        Initial water content assumed equal to field capacity
        """
        
        self.D = p['depth']
        self.Ksat = p['Ksat']
        self.pF = p['pF']
        
        self.Fc = psi_theta(self.pF, x=-1.0)
        self.Wp = psi_theta(self.pF, x=-150.0)
        
        self.poros = self.pF['ThetaS'] # porosity
        
        # water storages
        self.MaxSto = self.poros * self.D
        self.MaxPond = p['MaxPond']
        
        # initial state
        self.SurfSto = 0.0
        self.Wliq = self.Fc
        self.WatSto = self.D * self.Wliq
        self.Wair = self.poros - self.Wliq        
        self.Sat = self.Wliq / self.poros
        self.h = theta_psi(self.pF, self.Wliq) # water potential [m]
        self.Kh = self.Ksat * hydraulic_conductivity(self.pF, self.h)
        
        # relatively extractable water
        self.Rew = np.minimum( (self.Wliq - self.Wp) / (self.Fc - self.Wp + EPS), 1.0)
                          
    def run(self, dt, rr=0, et=0, latflow=0):
        """ Bucket model water balance
        Args: 
            dt [s]
            rr = potential infiltration [mm s-1 = kg m-2 s-1]
            et [mm s-1]
            latflow [mm s-1]
        Returns: 
            infil [mm] - infiltration [m]
            Roff [mm] - surface runoff
            drain [mm] - percolation /leakage
            roff [mm] - surface runoff
            et [mm] - evapotranspiration
            mbe [mm] - mass balance error
        """
        
        # fluxes        
        rr = rr * dt / WATER_DENSITY
        latflow = latflow * dt / WATER_DENSITY
        Qin = (rr + latflow) # m, potential inputs   
        et = et * dt / WATER_DENSITY
        
        # free drainage from profile
        drain = min(self.Kh * dt, max(0, (self.Wliq - self.Fc)) * self.D) # m
                
        # infiltr is restricted by availability, Ksat or available pore space
        infil = min(Qin, self.Ksat*dt, (self.MaxSto - self.WatSto + drain + et)) 

        # change in surface and soil water store
        SurfSto0 = self.SurfSto       
        dSto = (infil - drain - et)
        
        #in case of Qin excess, update SurfSto and route Roff
        q = Qin - infil
        if q > 0:
            dSur = min(self.MaxPond - self.SurfSto, q)
            roff = q - dSur #runoff, m
        else:
            dSur = 0.0
            roff = 0.0
        
        #update state variables
        self.update_state(dSto, dSur)
        
        #mass balance error
        mbe = dSto + (self.SurfSto - SurfSto0) - (rr + latflow - et - drain - roff)
        
        print(self.Rew)
        return WATER_DENSITY * infil, WATER_DENSITY * roff, WATER_DENSITY * drain, WATER_DENSITY * mbe 

    def update_state(self, dWat, dSur=0):
        """
        updates bucket state by computed dSto [m] and dSur [m] 
        """
        self.WatSto += dWat
        self.SurfSto +=dSur
        
        self.Wliq = self.poros * self.WatSto / self.MaxSto
        self.Wair = self.poros - self.Wliq
        self.Sat = self.Wliq/self.poros
        self.h = theta_psi(self.pF, self.Wliq)
        self.Kh = self.Ksat * hydraulic_conductivity(self.pF, self.h)            
        self.Rew = np.minimum((self.Wliq - self.Wp) / (self.Fc - self.Wp + EPS), 1.0)      
 

# soil water retention and hydraulic conductivity functions

def theta_psi(pF, x):
    """
    converts volumetric moisture content (m3m-3) to soil water potential
    Args:
        pF - vanGenuchten -model parameters (dict of floats)
        x - vol. water content (m3m-3, float)
    Returns:
        Psi (m)
    """
    # converts water content (m3m-3) to potential (m)
    
    ts = pF['ThetaS']
    tr = pF['ThetaR']
    a = pF['alpha']
    n = pF['n']
    m = 1.0 - np.divide(1.0, n)
    
    x = np.minimum(x, ts)
    x = np.maximum(x, tr)  # checks limits
    s = (ts - tr) / ((x - tr) + EPS)
    Psi = -1e-2 / a * (s**(1.0 / m) - 1.0)**(1.0 / n)  # m

    return Psi

def psi_theta(pF, x):
    """
    converts water potential (m) to water content (m3m-3)
        Args:
        pF - vanGenuchten -model parameters (dict of floats)
        x - soil water potential (m, float)
    Returns:
        Th - vol. water content (m3m-3, float)
    """
    x = 100 * np.minimum(x, 0)  # cm
    ts = pF['ThetaS']
    tr = pF['ThetaR']
    a = pF['alpha']
    n = pF['n']
    m = 1.0 - np.divide(1.0, n)
    
    Th = tr + (ts - tr) / (1.0 + abs(a * x)**n)**m
    return Th

def hydraulic_conductivity(pF, x, Ksat=1.0):
    # Hydraulic conductivity (vanGenuchten-Mualem)
    # x = water potential [m]
    x = 100 * np.minimum(x, 0)  # cm
    a = pF['alpha']
    n = pF['n']
    m = 1.0 - np.divide(1.0, n)    

    def relcond(x):
        Seff = 1.0 / (1.0 + abs(a * x)**n)**m
        r = Seff**0.5 * (1.0 - (1.0 - Seff**(1/m))**m)**2.0
        return r
    
    return Ksat * relcond(x)

# Evaporation algorithms

def eq_evap(AE, T, P=101300.0):
    """
    Calculates the equilibrium evaporation according to McNaughton & Spriggs,\
    1986.
    INPUT:
        AE - Available energy (Wm-2)
        T - air temperature (degC)
        P - pressure (Pa)
    OUTPUT:
        equilibrium evaporation rate (Wm-2)
        lat. heat of vaporization (J kg-1)
    """
    NT = 273.15
    cp = 1004.67  # J/kg/K
    
    # L lat. heat of vaporization (J/kg), esa = saturation vapor pressure (Pa), 
    # slope of esa (Pa K-1), psychrometric constant g (Pa K-1)
    L = 1e3 * (3147.5 - 2.37 * (T + NT))  # lat heat of vapor [J/kg]
    esa = 1e3 * (0.6112 * np.exp((17.67 * T) / (T + 273.16 - 29.66)))  # Pa

    s = 17.502 * 240.97 * esa / ((240.97 + T) ** 2)
    g = P * cp / (0.622 * L)
    
    # equilibrium evaporation # Wm-2 = Js-1m-2
    x = np.divide((AE * s), (s + g))  

    x = np.maximum(x, 0.0)
    return x, L

# def CanopyConductance(p, LAI, Qp, Ta, D, fPheno=1.0, Rew=1.0, CO2=380.0):
#     """
#     Canopy conductance, follows Launiainen et al. 2019 HESS
#     Args:
#         p - parameters (dict)
#         LAI - 1-sided leaf-area index (m2m-2)
#         Qp - incoming PAR (W m-2)
#         Ta - air temperature (degC)
#         D - vapor pressure deficit (kPa)
#         fPheno - seasonal cycle modifier (-)
#         Rew - relatively extractable soil water (-)
#         CO2 - ambient CO2 (ppm)
#     Returns:
#         Gc - canopy conductance (canopy-integrated stomatal conductance)  (m s-1)
#     """
#     D = np.max(D, EPS) # VPD alwas >= 0
    
#     rhoa = 101300.0 / (8.31 * (Ta + 273.15)) # mol m-3
    
#     Amax = p['Amax'] # light-saturated leaf photosynthetic rate (umol m-2 s-1)
#     g1 = p['g1'] # USO model parameter (unit?)
#     kp = p['kp']  # (-) attenuation coefficient for PAR
#     q50 = p['q50']  # Wm-2, half-sat. of leaf light PAR response
#     rw = p['rw']  # rew parameter; start of conductance decay
#     rwmin = p['rwmin']  # rew parameter; minimum relative conductance

#     # leaf level light-saturated stomatal conductance gs (m/s)
#     gs = 1.6*(1.0 + g1 / np.sqrt(D)) * Amax / CO2 / rhoa

#     """--- environmental modifiers """

#     # fQ (canopy-level light response): Saugier & Katerji, 1991 Agric. For. Met., eq. 4. Leaf light response = Qp / (Qp + q50)
#     fQ = 1./ kp * np.log((kp*Qp + q50) / (kp*Qp*np.exp(-kp * LAI) + q50 + EPS) )

#     # the next formulation is from Leuning et al., 2008 WRR for daily Gc; they refer to 
#     # Kelliher et al. 1995 AFM but the resulting equation is not exact integral of K95.        
#     # fQ = 1./ kp * np.log((Qp + q50) / (Qp*np.exp(-kp*self.LAI) + q50))

#     # soil moisture response
#     fRew = np.minimum(1.0, np.maximum(Rew / rw, rwmin))

#     # CO2 -response of canopy conductance, derived from APES-simulations
#     # (Launiainen et al. 2016, Global Change Biology). relative to 380 ppm
#     fCO2 = 1.0 - 0.387 * np.log(CO2 / 380.0)
    
#     # canopy conductance
#     Gc = gs * fQ * fRew * fCO2 * fPheno
    
#     return Gc
    
