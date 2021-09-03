import numpy as np
from atmosredox import *
from numba import njit, jit, float32

def GCO(T):
    dG = 0.05829 * T - 252.033
    return dG * 1e3

def CO_Cratio(fO2, T):
    muCO = GCO(T)
    muC = 0
    muO2 = 0

    dG = muC + 0.5 * muO2 - muCO
    Keq = np.exp(dG / (R * T))

    xCO_xC = Keq * np.sqrt(fO2)
    return xCO_xC

def GCO2(T):
    dG = 7.43e-4 * T - 397.691
    return dG * 1e3

def CO2_COratio(fO2, T):
    muCO2 = GCO2(T)
    muCO = GCO(T)
    muO2 = 0

    dG = muCO + 0.5 * muO2 - muCO2
    Keq = np.exp(dG / (R * T))
    XCdioxide_CO = Keq * np.sqrt(fO2)
    return XCdioxide_CO

def GCH4(T):
    dG = 0.1115725 * T - 92.3855
    return dG * 1e3

def COH2O_CH4ratio(fO2, T):
    muCO = GCO2(T)
    muCH4 = GCH4(T)
    muH2O = GH2O(T)

    dG = muCH4 - muCO - 2 * muH2O
    Keq = np.exp(dG / (R * T))
    Xspecies_XCH4 = Keq * fO2**(3 / 2)
    return Xspecies_XCH4

@jit(nopython=True)
def GFeO(T):
    dG = 0.05277 * T - 254.1475
    return dG

@jit(nopython=True)
def Keq_FeO_H2O(T):
    Keq = np.exp((GH2O(T) - GFeO(T)) / R / T)
    return Keq

@jit(nopython=True)
def Keq_FeO_CO2(T):
    Keq = np.exp((GCO2(T) - GFeO(T)) / R / T)
    return Keq
