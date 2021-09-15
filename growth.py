from equilibrium import *
import pandas as pd
import numpy as np

# knowns
r_planet = 100 * 1e3
k = 3.926763924239811e-12 # using Mars' mass-g relationship
G = 6.674e-11


# constants for metal pond timescale equation
rho_2 = 8000 # kg m^-3
mu_1 = 3 * 1e19 # Pa s
mu_2 = 10 # Pa s
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


def calculate_g(M_p):
    """Calculate gravitational acceleration using g ~ M_p^(0.503) and M_Mars = 6.39e23 kg, g_Mars = 3.7 m / s"""
    return M_p**0.503 * k

def calculate_h(melt_vol, r_planet):
    """Returns the depth of mantle that is melted by an impactor."""
    return r_planet - (r_planet**3 - 3 / 4 / math.pi * melt_vol)**(1 / 3)

def Peq(g, h):
    """Returns P_eq in GPa."""
    return rho_mantle * g * h / 1e9

def mantle_density(Peq):
    """Already written"""
    if Peq < 24e9:
        rho_0 = 3300
        K0d = 5
        K0 = 125e9
    
    else:
        rho_0 = 4100
        K0d = 3.97
        K0 = 247e9

    return rho_0 * (1 + K0d * Peq / K0)**(1 / K0d)

def Peq2(h, r_planet, m_planet):
    """Already written"""
    N = 100
    z = np.zeros(N)
    r = np.zeros(N)
    p = np.zeros(N)
    rho = np.zeros(N)
    g = np.zeros(N)

    dr = h / (N - 1)
    r[0] = r_planet
    for i in range(1, N):
        z[i] = h * i / (N -1)
        r[i] = r_planet - z[i]
        rho[i] = mantle_density(p[i - 1])
        m_planet -= rho[i] * (4 * math.pi * r[i]**2 * dr)

        # using -dP/dr = rho * m * g
        g[i] = G * m_planet / r[i]**2
        p[i] = p[i - 1] + rho[i] * g[i] * dr

    return p[-1]

def added_depth(initial_radius, added_mass, rho):
    """Return increase in radius / depth of core / mantle respectively given added mass."""
    return (3 / 4 / math.pi * added_mass / rho + initial_radius**3)**(1 / 3) - initial_radius

def core_radius(body_radius):
    """Returns initial core radius upon the assumption that initial metal mass is 25% of the total mass."""
    return (body_radius**3 / 9)**(1 / 3)

def core_radius2(body_radius):
    """Returns initial core radius upon the assumption that initial metal mass is 20% of the total mass."""
    return (body_radius**3 / 9.4)**(1 / 3)

def body_mass(r_body, core_radius):
    return 4 / 3 * math.pi * (core_radius**3 * rho_core + (r_body - core_radius)**3 * rho_mantle)

def ocean_mass(r_planet, depth):
    return 4 / 3 * math.pi * (r_planet**3 - (r_planet - depth)**3) * rho_mantle

def convert_moles_to_mass(mol_element, element_molar_mass):
    """Returns the mass in kg of an element given its molar mass and number of moles."""
    return mol_element * element_molar_mass / 1000

def Teq(Peq):
    """Input Peq in Pa."""
    return 0.4 * 1661.2 * (Peq / 1.336 / 10**9 + 1)**(1 / 7.437) + 0.6 * 1982.1 * (Peq / 6.594 / 1e9 + 1)**(1 / 5.374)

def Teq2(Peq):
    """Input Peq in GPa."""
    if Peq < 24:
        return 1874 + 55.43 * Peq - 1.74 * Peq**2 + 0.0193 * Peq**3
    else:
        return 1249 + 58.28 * Peq - 0.395 * Peq**2 + 0.0011 * Peq**3

def sphere_radius(mass, density):
    return (3 / 4/ math.pi * mass / density)**(1 / 3)

def shell_width(mass, density, inner_radius):
    return (3 / 4 / math.pi * mass / density + inner_radius**3)**(1 / 3) - inner_radius

def shell_width_from_outer(mass, density, outer_radius):
    return outer_radius - (- 3 / 4 / math.pi * mass / density + outer_radius**3)**(1 / 3)

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