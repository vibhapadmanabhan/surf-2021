from numba.core.errors import deprecated
from equilibrium import *
import pandas as pd
from planetesimal_initial_parameters import *

# lists

# chemical
# core
X_Fe = []
X_Si = []
X_Ni = []
X_Va = []
# mantle
X_FeO = []
X_SiO2 = []
X_NiO = []
X_Mg = []
X_VO = []
X_FeO_impactor = []
# physical 
gravity = []
pressure = []
temperature = []
planet_size = []
impactor_size = []
mantle_depth = []
fO2 = []
metal_pond_lifetime = []


def save_data(X_Fe, X_Si, X_Ni, X_Va, X_FeO, X_SiO2, X_NiO, X_Mg, X_VO, X_FeO_impactor, gravity, pressure, temperature, planet_size, impactor_size, mantle_depth, fO2, filename):
    df = pd.DataFrame()
    df["X_Fe"] = X_Fe
    df["X_Si"] = X_Si
    df["X_Ni"] = X_Ni
    df["X_V"] = X_Va
    df["X_FeO"] = X_FeO
    df["X_SiO2"] = X_SiO2
    df["X_NiO"] = X_NiO
    df["X_Mg"] = X_Mg
    df["X_V2O3"] = X_VO
    df["g (m/s)"] = gravity
    df["P (GPa)"] = pressure
    df["T (K)"] = temperature
    df["Planet radius (km)"] = planet_size
    df["Impactor size (km)"] = impactor_size
    df["Magma ocean depth (km)"] = mantle_depth
    df["ln(fO2)_IW"] = fO2
    df["FeO wt\% in impactor"] = X_FeO_impactor
    df.to_csv(filename, sep='\t', index=False)

def metal_pond_timescale(r1, r2, mantle_solid_depth, planet_core_radius):
    """Return metal pond timescale in kyr."""
    solid_mantle_core_mass = calculate_vol(planet_core_radius) * 8000 + calculate_shell_vol(mantle_solid_depth, planet_core_radius + mantle_solid_depth) * 3000
    rho_1 = solid_mantle_core_mass / calculate_vol(planet_core_radius + mantle_solid_depth)
    sigma = 10**(-10.13) * r1**1.86 * (rho_2 - rho_1)**1.3 * rho_1**0.87 * (r2 - r1)**0.047 * mu_1**(-0.98) * mu_2**(-0.054)
    return 1 / sigma / 3600 / 24 / 365 / 1000