import numpy as np
from numba import njit, jit, float32

R   = 8.31446  # universal gas constant


def fO2_fromIW(fO2_dIW, T):
    ''' convert fO2 in terms of [bar]. The relation ln(fO2) = 2 ln(xFeO/xFe) shows the value in terms of [difference from IW] (by Yoshi Miyazaki)'''

    # log10(fO2) is -12 at 1200 degC, -8 at 1700 degC.
    # assume a log-linear relationship for now
    fO2_IW = 10**(-12 + (T-1473.15)/500*4.)

    # convert
    fO2    = fO2_IW * fO2_dIW
    return fO2


# Gibbs free energies of formation fit with a linear regression using data on JANAF database
@jit(nopython=True)
def GCO(T):
    dG = 0.05829 * T - 252.033
    return dG * 1e3


@jit(nopython=True)
def GCO2(T):
    dG = 7.43e-4 * T - 397.691
    return dG * 1e3


@jit(nopython=True)
def GH2O(Tin):
    ''' dG (J/mol) as a function of T. Data taken from JANAF database and fit using a linear regression (by Yoshi Miyazaki)'''
    dG = 0.0547*Tin - 246.56   # in kJ/mol. Tin is in [K].

    return dG*1e3  # in J/mol.


@jit(nopython=True)
def GCH4(T):
    dG = 0.1115725 * T - 92.3855
    return dG * 1e3


@jit(nopython=True)
def GFeO(T):
    dG = 0.05277 * T - 254.1475
    return dG * 1e3


# oxidized : reduced compound ratios calculated using fugacity
def CO_Cratio(fO2, T):
    # C + 0.5 O2 <--> CO
    muCO = GCO(T)
    muC = 0
    muO2 = 0

    dG = muCO - muC - 0.5 * muO2
    Keq = np.exp(-dG / (R * T))

    xCO_xC = Keq * np.sqrt(fO2)
    return xCO_xC


def CO2_COratio(fO2, T):
    # CO + 0.5 O2 <--> CO2
    muCO2 = GCO2(T)
    muCO = GCO(T)
    muO2 = 0

    dG = muCO2 - muCO - 0.5 * muO2
    Keq = np.exp(-dG / (R * T))
    XCdioxide_CO = Keq * np.sqrt(fO2)
    return XCdioxide_CO


def CO2H2O_CH4ratio(fO2, T):
    # """Calculates equilibrium constant for the reaction CH4 + 2 O2 <--> CO2 + 2 H2O"""
    muCO2 = GCO2(T)
    muCH4 = GCH4(T)
    muH2O = GH2O(T)

    dG = muCO2 + 2 * muH2O - muCH4
    Keq = np.exp(-dG / (R * T))
    Xspecies_XCH4 = Keq * fO2**2
    return Xspecies_XCH4


def H2O_H2ratio(fO2, Tin):
    ''' solve for H2O/H2 ratio using equilibrium constant (by Yoshi Miyazaki)'''
    muH2O = GH2O(Tin)
    muO2  = 0.
    muH2  = 0.
    
    dG  = - muH2 - 0.5*muO2 + muH2O  # in J/mol
    Keq = np.exp(-dG/(R*Tin))       # convert to equilibrium constant
    xH2O_xH2 = Keq*np.sqrt(fO2)    # equilibrium relation
    
    return xH2O_xH2

def Keq_H2O_H2(Tin):
    muH2O = GH2O(Tin)
    muO2  = 0.
    muH2  = 0.
    dG  = - muH2 - 0.5*muO2 + muH2O  # in J/mol
    Keq = np.exp(-dG/(R*Tin))  

    return Keq


# equilibrium constants for re-equilibrium
@jit(nopython=True)
def Keq_FeO_H2O(T):
    """Calculates equilibrium constant for the reaction FeO + H2 <--> Fe + H2O"""
    Keq = np.exp(-(GH2O(T) - GFeO(T)) / R / T)
    return Keq

@jit(nopython=True)
def Keq_FeO_CO2(T):
    """Calculates equilibrium constant for the reaction FeO + CO <--> Fe + CO2"""
    Keq = np.exp(-(GCO2(T) - GFeO(T) - GCO(T)) / R / T)
    return Keq

@jit(nopython=True)
def Keq_FeO_CH4(T):
    """Calculates equilibrium constant for the reaction FeO + CH4 <--> Fe + CO + 2 H2"""
    Keq = np.exp(-(GCO(T) - GFeO(T) - GCH4(T)) / R / T)
    return Keq
