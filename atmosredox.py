# atmosredox.py
# Aug 12, 2021 (Yoshi Miyazaki)

import numpy             as np
import matplotlib.pyplot as plt
from numba import jit, njit, float32

# import matplotlib.patches as ptch
# from   matplotlib import rc
# from   matplotlib.font_manager import FontProperties
# from   matplotlib.ticker       import MultipleLocator, FormatStrFormatter
# import sys

# plt.rcParams['figure.figsize']= 4*1.414, 4*1
# plt.rcParams['font.family']='sans-serif'
# plt.rcParams['font.size']     = 14
# plt.rcParams['axes.linewidth']= 1.0
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text',usetex=True)
# rc('text.latex', preamble=r'\usepackage{sfmath}')

R   = 8.31446  # universal gas constant

def fO2_fromIW(fO2_dIW, T):
    ''' convert fO2 in terms of [bar]. The relation ln(fO2) = 2 ln(xFeO/xFe) shows the value in terms of [difference from IW] '''

    # log10(fO2) is -12 at 1200 degC, -8 at 1700 degC.
    # assume a log-linear relationship for now
    fO2_IW = 10**(-12 + (T-1473.15)/500*4.)

    # convert
    fO2    = fO2_IW * fO2_dIW

    return fO2

@jit(nopython=True)
def GH2O(Tin):
    ''' dG (J/mol) as a function of T. Data taken from JANAF database and fit using a linear regression '''
    dG = 0.0547*Tin - 246.56   # in kJ/mol. Tin is in [K].

    return dG*1e3  # in J/mol.

def H2O_H2ratio(fO2, Tin):
    ''' solve for H2O/H2 ratio using equilibrium constant '''
    muH2O = GH2O(Tin)
    muO2  = 0.
    muH2  = 0.
    
    dG  = - muH2 - 0.5*muO2 + muH2O  # in J/mol
    Keq = np.exp(-dG/(R*Tin))       # convert to equilibrium constant
    xH2O_xH2 = Keq*np.sqrt(fO2)    # equilibrium relation
    
    return xH2O_xH2

# calculate H2O/H2 ratio for a range of oxygen fugacity
# set temperature
T = 1800

# prepare array of fO2 and H2/(H2+H2O) ratio
l_fO2 = np.array([-3., -2, -1, 0, 1, 2, 3, 4])
l_rat = np.zeros(len(l_fO2))

for i, fO2_dIW in enumerate(l_fO2):
    
    fO2 = fO2_fromIW(np.exp(fO2_dIW), T)  # convert ln(fO2(dIW)) into fO2 [Pa]
    r   = H2O_H2ratio(fO2, T)
    
    l_rat[i] = 1/(1+r)


# create a plot
fig, ax = plt.subplots()
ax.plot(l_fO2, l_rat, color='k')
ax.set_xlim([np.min(l_fO2), np.max(l_fO2)])
ax.set_ylim([0, 1])
ax.set(xlabel="Oxygen fugacity [ln(fO2/IW)", ylabel="H$_2$/(H$_2$+H$_2$O) at 1800 K, 1 bar")

# plot options
ax.tick_params(direction='in', length=3,   which='major', top=True, right=True)
ax.tick_params(direction='in', length=1.5, which='minor', top=True, right=True)
plt.tight_layout()
plt.savefig("./fig_fO2_H2ratio.pdf", transparent=True)
