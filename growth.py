from equilibrium import *
import pandas as pd

# knowns
r_planet = 100 * 1e3
k = 3.926763924239811e-12 # using Mars' mass-g relationship

# starting mantle composition of oxygen-free elements, and oxygen, in wt% (assuming FeO is 8 wt%, MgO is 42 wt%, SiO2 is 50 - 0.00606 wt%)
fe_s = 6.2172
mg_s = 25.325
si_s = 23.33
ni_s = 0
v_s = 0.004119
o_s = 100 - fe_s - mg_s - si_s - ni_s - v_s

# starting core composition in wt%
fe = 85
ni = 15
si = 0
v = 0

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

# physical 
gravity = []
pressure = []
temperature = []
planet_size = []
impactor_size = []
mantle_depth = []
fO2 = []


def calculate_g(M_p):
    """Calculate gravitational acceleration using g ~ M_p^(0.503) and M_Mars = 6.39e23 kg, g_Mars = 3.7 m / s"""
    return M_p**0.503 * k

def calculate_h(melt_vol, r_planet):
    """Returns the depth of mantle that is melted by an impactor."""
    return r_planet - (r_planet**3 - 3 / 4 / math.pi * melt_vol)**(1 / 3)

def Peq(g, h):
    """Returns P_eq in GPa."""
    return rho_mantle * g * h / 1e9

def added_depth(initial_radius, added_mass, rho):
    """Return increase in radius / depth of core / mantle respectively given added mass."""
    return (3 / 4 / math.pi * added_mass / rho + initial_radius**3)**(1 / 3) - initial_radius

def core_radius(body_radius):
    """Returns initial core radius upon the assumption that initial metal mass is 25% of the total mass."""
    return (body_radius**3 / 9)**(1 / 3)

def body_mass(r_body, core_radius):
    return 4 / 3 * math.pi * (core_radius**3 * rho_core + (r_body - core_radius)**3 * rho_mantle)

def ocean_mass(r_planet, depth):
    return 4 / 3 * math.pi * (r_planet**3 - (r_planet - depth)**3) * rho_mantle

def convert_moles_to_mass(mol_element, element_molar_mass):
    """Returns the mass in kg of an element given its molar mass and number of moles."""
    return mol_element * element_molar_mass / 1000

def Teq(Peq):
    return 0.4 * 1661.2 * (Peq / 1.336 / 10**9 + 1)**(1 / 7.437) + 0.6 * 1982.1 * (Peq / 6.594 / 1e9 + 1)**(1 / 5.374)

def sphere_radius(mass, density):
    return (3 / 4/ math.pi * mass / density)**(1 / 3)

def shell_width(mass, density, inner_radius):
    return (3 / 4 / math.pi * mass / density + inner_radius**3)**(1 / 3) - inner_radius

def save_data(X_Fe, X_Si, X_Ni, X_Va, X_FeO, X_SiO2, X_NiO, X_VO, gravity, pressure, temperature, planet_size, impactor_size, mantle_depth, fO2, filename):
    df = pd.DataFrame()
    df["X_Fe"] = X_Fe
    df["X_Si"] = X_Si
    df["X_Ni"] = X_Ni
    df["X_V"] = X_Va
    df["X_FeO"] = X_FeO
    df["X_SiO2"] = X_SiO2
    df["X_NiO"] = X_NiO
    df["X_V2O3"] = X_VO
    df["g (m/s)"] = gravity
    df["P (GPa)"] = pressure
    df["T (K)"] = temperature
    df["Planet radius (km)"] = planet_size
    df["Impactor size (km)"] = impactor_size
    df["Magma ocean depth (km)"] = mantle_depth
    df["ln(fO2)_IW"] = fO2
    df.to_csv(filename, sep='\t', index=False)

    # with open("./data/rapid-solidification/X_FeO.txt", "w") as f:
    #     for val in X_FeO:
    #         f.write("%s\n" % val)

    # with open("./data/rapid-solidification/X_Fe.txt", "w") as f:
    #     for val in X_Fe:
    #         f.write("%s\n" % val)

    # with open("./data/rapid-solidification/X_Va.txt", "w") as f:
    #     for val in X_Va:
    #         f.write("%s\n" % val)

    # with open("./data/rapid-solidification/X_V2O3.txt", "w") as f:
    #     for val in X_V2O3:
    #         f.write("%s\n" % val)

    # with open("./data/rapid-solidificationX_Ni.txt", "w") as f:
    #     for val in X_Ni:
    #         f.write("%s\n" % val)

    # with open("./data/rapid-solidification/X_NiO.txt", "w") as f:
    #     for val in X_NiO:
    #         f.write("%s\n" % val)

    # with open("./data/rapid-solidification/X_SiO2.txt", "w") as f:
    #     for val in X_SiO2:
    #         f.write("%s\n" % val)

    # with open("./data/rapid-solidification/X_Si.txt", "w") as f:
    #     for val in X_Si:
    #         f.write("%s\n" % val)

    # with open("./data/rapid-solidification/temperature.txt", "w") as f:
    #     for val in temperature:
    #         f.write("%s\n" % val)

    # with open("./data/rapid-solidification/pressure.txt", "w") as f:
    #     for val in pressure:
    #         f.write("%s\n" % val)

    # with open("./data/rapid-solidification/gravity.txt", "w") as f:
    #     for val in gravity:
    #         f.write("%s\n" % val)

    # with open("./data/rapid-solidification/h.txt", "w") as f:
    #     for val in mantle_depth:
    #         f.write("%s\n" % val)
    
    # with open("./data/rapid-solidification/planet_size.txt", "w") as f:
    #     for val in planet_size:
    #         f.write("%s\n" % val)

    # with open("./data/rapid-solidification/fO2.txt", "w") as f:
    #     for val in fO2:
    #         f.write("%s\n" % val)


def save_deep_MO_data(X_FeO, X_Fe, X_SiO2, X_Si, X_Va, X_V2O3, X_Ni, X_NiO, pressure, temperature, gravity, planet_size, fO2, mantle_depth):
    with open("./data/deep-MO/deep_MO_X_FeO.txt", "w") as f:
        for val in X_FeO:
            f.write("%s\n" % val)

    with open("./data/deep-MO/deep_MO_X_Fe.txt", "w") as f:
        for val in X_Fe:
            f.write("%s\n" % val)

    with open("./data/deep-MO/deep_MO_X_Va.txt", "w") as f:
        for val in X_Va:
            f.write("%s\n" % val)

    with open("./data/deep-MO/deep_MO_X_V2O3.txt", "w") as f:
        for val in X_V2O3:
            f.write("%s\n" % val)

    with open("./data/deep-MO/deep_MO_X_Ni.txt", "w") as f:
        for val in X_Ni:
            f.write("%s\n" % val)

    with open("./data/deep-MO/deep_MO_X_NiO.txt", "w") as f:
        for val in X_NiO:
            f.write("%s\n" % val)

    with open("./data/deep-MO/deep_MO_X_SiO2.txt", "w") as f:
        for val in X_SiO2:
            f.write("%s\n" % val)

    with open("./data/deep-MO/deep_MO_X_Si.txt", "w") as f:
        for val in X_Si:
            f.write("%s\n" % val)

    with open("./data/deep-MO/deep_MO_temperature.txt", "w") as f:
        for val in temperature:
            f.write("%s\n" % val)

    with open("./data/deep-MO/deep_MO_pressure.txt", "w") as f:
        for val in pressure:
            f.write("%s\n" % val)

    with open("./data/deep-MO/deep_MO_gravity.txt", "w") as f:
        for val in gravity:
            f.write("%s\n" % val)

    with open("./data/deep-MO/deep_MO_h.txt", "w") as f:
        for val in mantle_depth:
            f.write("%s\n" % val)
    
    with open("./data/deep-MO/deep_MO_planet_size.txt", "w") as f:
        for val in planet_size:
            f.write("%s\n" % val)

    with open("./data/deep-MO/deep_MO_fO2.txt", "w") as f:
        for val in fO2:
            f.write("%s\n" % val)

    