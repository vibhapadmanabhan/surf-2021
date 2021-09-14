import numpy as np
from atmosredox import *
from numba import njit, jit, float32

def Keq_H2_H2O(T):
    Keq = np.exp(-GH2O / R / T)
    return Keq

def Keq_CH4(T):
    Keq = np.exp(-(2 * GH2O - GCO2 - GCH4) / R / T)
    return Keq

def Keq_O2_CO2(T):
    Keq = np.exp(-GCO2 / R / T)
    return Keq
    
@jit(nopython=True)
def GCO(T):
    dG = 0.05829 * T - 252.033
    return dG * 1e3

def CO_Cratio(fO2, T):
    # C + 0.5 O2 <--> CO
    muCO = GCO(T)
    muC = 0
    muO2 = 0

    dG = muCO - muC - 0.5 * muO2
    Keq = np.exp(-dG / (R * T))

    xCO_xC = Keq * np.sqrt(fO2)
    return xCO_xC

@jit(nopython=True)
def GCO2(T):
    dG = 7.43e-4 * T - 397.691
    return dG * 1e3

def CO2_COratio(fO2, T):
    # CO + 0.5 O2 <--> CO2
    muCO2 = GCO2(T)
    muCO = GCO(T)
    muO2 = 0

    dG = muCO2 - muCO - 0.5 * muO2
    Keq = np.exp(-dG / (R * T))
    XCdioxide_CO = Keq * np.sqrt(fO2)
    return XCdioxide_CO

@jit(nopython=True)
def GCH4(T):
    dG = 0.1115725 * T - 92.3855
    return dG * 1e3

def CO2H2O_CH4ratio(fO2, T):
    # CH4 + 2O2 <--> CO2 + 2H2O
    muCO2 = GCO2(T)
    muCH4 = GCH4(T)
    muH2O = GH2O(T)

    dG = muCO2 + 2 * muH2O - muCH4
    Keq = np.exp(-dG / (R * T))
    Xspecies_XCH4 = Keq * fO2**2
    return Xspecies_XCH4

@jit(nopython=True)
def GFeO(T):
    dG = 0.05277 * T - 254.1475
    return dG

@jit(nopython=True)
def Keq_FeO_H2O(T):
    # FeO + H2 <--> Fe + H2O
    Keq = np.exp(-(GH2O(T) - GFeO(T)) / R / T)
    return Keq

@jit(nopython=True)
def Keq_FeO_CO2(T):
    # FeO + CO <--> Fe + CO2
    Keq = np.exp(-(GCO2(T) - GFeO(T) - GCO(T)) / R / T)
    return Keq

@jit(nopython=True)
def Keq_FeO_CH4(T):
    # FeO + CH4 <--> Fe + CO + 2H2
    Keq = np.exp(-(GCO(T) - GFeO(T) - GCH4(T)) / R / T)
    return Keq


