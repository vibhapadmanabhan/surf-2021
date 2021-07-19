import math
import numpy as np 
import matplotlib.pyplot as plt

# constants
r_planet = 1000 * 1E3 # m
rho_mantle = 3000
rho_core = 8000
P_eq = 0.7 * 25 # in GPa, assuming P_eq = 0.7 * P_cmb. CMB pressure depends on size of the planet.
T_cmb = 4000

# Mantle composition: chondritic, in wt%, ignoring trace elements and Ni, Fe.
si = 10.7
fe_s = 50
v = 0.00606

# Core composition, in wt%
fe = 85
ni = 15

# molar masses (g / mol)
molar_mass_fe = 55.8
molar_mass_ni = 58.7
molar_mass_si = 28
molar_mass_o = 16
molar_mass_v = 50.9415

planet_mantle_depth = 300 * 1E3  # assumed

def calculate_impactor_core_radius(impactor_radius):
    """
    Returns the radius of the core of a body in meters using Earth's core-mantle mass ratio.
    """
    return (0.15 * impactor_radius**3)**(1 / 3)


def calculate_shell_vol(mantle_depth, r_obj):
    """Returns the volume of a spherical shell given the width of the shell and its outer radius."""
    return 4 / 3 * math.pi * (r_obj**3 - (r_obj - mantle_depth)**3)


def calculate_vol(r_obj):
    """Returns the volume of a sphere."""
    return 4 / 3 * math.pi * r_obj**3


def calculate_element_mass(section_mass, element_wt_percent):
    """Returns the mass of an element in g in the mantle/core given the wt% of the element. """
    return section_mass * element_wt_percent * 0.01 * 1000


def convert_mass_to_moles(compound_mass, compound_molar_mass):
    """Converts mass to moles given the mass of a compound in g and its molar mass in g / mol."""
    return compound_mass / compound_molar_mass


def calculate_total_element_moles(element, element_molar_mass, planet_mantle_depth, r_planet, r_impactor, impactor_core_radius):
    """Given an impactor and planet, calculates the total number of moles of Si, Ni, Fe, or FeO.
    
    Does not calculate moles of O."""
    if element == si or element == fe_s:
        mass = (calculate_shell_vol(planet_mantle_depth, r_planet) + calculate_shell_vol(r_impactor - impactor_core_radius, r_impactor)) * rho_mantle
        element_mass = calculate_element_mass(mass, element)
        return convert_mass_to_moles(element_mass, element_molar_mass)
    mass = calculate_vol(impactor_core_radius) * rho_core
    element_mass = calculate_element_mass(mass, element)
    return convert_mass_to_moles(element_mass, element_molar_mass)


def calculate_kd(metal, T, P, a, b, c):
    """Calculates K_d using constants given in Rubie 2015 Table S1 and Rubie 2015 equation (3)."""
    if metal == "si":
        return math.exp(a + b / T + (-155 * P + 2.26 * P**2 - 0.011 * P**3) / T)
    else:
        return math.exp(a + b / T + c * P / T)


# actual distribution constants from Rubie
actual_kd_si = calculate_kd("si", T_cmb, P_eq, 2.98, -15934, None)
actual_kd_ni = calculate_kd("ni", T_cmb, P_eq, 1.06, 1553, -98)
actual_kd_v = calculate_kd("v", T_cmb, P_eq, -0.48, -5063, 0)

def f(fe_metal, mol_fe, mol_ni, mol_si, mol_o, mol_v, v_metal, P_eq):
    """
    Returns the difference between the actual value of K_d(Si-Fe) and calculated value for K_d(Si-Fe). (Root of function will be found at the 
    correct value of K_d)
    """
    actual_kd_si = calculate_kd("si", T_cmb, P_eq, 2.98, -15934, None)
    fe_sil = mol_fe - fe_metal
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil) / 2  # ignore mass of V oxides
    si_metal = mol_si - si_sil
    conc_si_metal = si_metal / (si_metal + fe_metal + ni_metal + v_metal)
    conc_si_sil = si_sil / (si_sil + fe_sil + ni_sil + (mol_v - v_metal))
    conc_fe_metal = fe_metal / (fe_metal + ni_metal + si_metal + v_metal)
    conc_fe_sil = fe_sil / (fe_sil + ni_sil + si_sil + (mol_v - v_metal))
    return  conc_si_metal * conc_fe_sil**2 / conc_si_sil / conc_fe_metal**2 - actual_kd_si


# still working on this function
def g(fe_metal, mol_fe, mol_ni, mol_si, mol_o, P_eq, ni_metal, si_metal, mol_v, v_metal):
    """
    Returns the difference between the actual value of K_d(V-Fe) and calculated value for K_d(V-Fe). (Root of function will be found at the 
    correct value of K_d)
    """
    actual_kd_v = calculate_kd("v", T_cmb, P_eq, -0.48, -5063, 0)
    fe_sil = mol_fe - fe_metal
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil) / 2  # ignore mass of V oxides
    si_metal = mol_si - si_sil
    conc_v_metal = v_metal / (si_metal + fe_metal + ni_metal + v_metal)
    conc_v_sil = (mol_v - v_metal) / (si_sil + fe_sil + ni_sil + mol_v - v_metal)
    conc_fe_metal = fe_metal / (fe_metal + ni_metal + si_metal + v_metal)
    conc_fe_sil = fe_sil / (fe_sil + ni_sil + si_sil + mol_v - v_metal)
    return  conc_v_metal * conc_fe_sil**(3 / 2) / conc_v_sil / conc_fe_metal**(3 / 2) - actual_kd_v


def root_bracket(func, mol_fe, mol_ni, mol_o, mol_si, P_eq, ni_metal, si_metal, mol_v, v_metal):
    """Finds the lower bound of an interval that can be used for the bisection search."""
    val = 10
    if func == f:
        FB = func(mol_fe, mol_fe, mol_ni, mol_si, mol_o, P_eq)
        while val <= mol_fe and FB * f(val, mol_fe, mol_ni, mol_si, mol_o, P_eq) > 0:
            val *= 1.2  # multiply by small enough number for root to be found
        return val
    elif func == g:
        FB = g(mol_fe, mol_fe, mol_ni, mol_si, mol_o, P_eq, ni_metal, si_metal, mol_v, v_metal)
        while val <= mol_fe and FB * g(val, mol_fe, mol_ni, mol_si, mol_o, P_eq, ni_metal, si_metal, mol_v, v_metal) > 0:
            val *= 1.2  # multiply by small enough number for root to be found
        return val


def bisection_search(func, a, b, eps, mol_fe, mol_ni, mol_si, mol_o, P_eq, ni_metal, si_metal, mol_v, v_metal):
    while True:
        FA = func(a, mol_fe, mol_ni, mol_si, mol_o, P_eq, ni_metal, si_metal, mol_v, v_metal)
        elem_metal = (a + b) / 2
        FP = func(fe_metal, mol_fe, mol_ni, mol_si, mol_o, P_eq, ni_metal, si_metal, mol_v, v_metal)
        if np.abs(FP) <= eps: # close enough to true value
            break
        if FA * FP > 0:
            a = elem_metal
        else:
            b = elem_metal
    return elem_metal


def calculate_ln_o_iw_fugacity(X_FeO, X_Fe):
     return 2 * math.log(X_FeO / X_Fe)


# still working on everything below, feel free to comment out
# r_impactor = [i * 10 * 1e3 for i in range(1, 10)]
# # r_impactor = 100 * 1e3
# # P_eq = [i * 0.1 * P_cmb for i in range(1, 10)]
# X_FeO = []
# X_Fe = []
# X_Si = []
# X_V2O3 = []
# v_metal = 0
# # mantle_depth = [i * 1e3 for i in range(300, 701)]
# for i in range(len(r_impactor)):
#     mol_si = calculate_total_element_moles(si, molar_mass_si, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
#     mol_fe_s = calculate_total_element_moles(fe_s, molar_mass_fe + molar_mass_o, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
#     mol_fe = calculate_total_element_moles(fe, molar_mass_fe, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i])) + mol_fe_s
#     mol_ni = calculate_total_element_moles(ni, molar_mass_ni, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
#     # consider vanadium in the mantle of the planet
#     mol_v = calculate_total_element_moles(v, molar_mass_v, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))

#     mol_o = 2 * mol_si + mol_fe_s

#     fe_metal = bisection_search(f, root_bracket(f, mol_fe, mol_ni, mol_si, mol_o, P_eq, None, None, None, None), mol_fe, 10e-7, mol_fe, mol_ni, mol_si, mol_o, P_eq)
    
#     fe_sil = mol_fe - fe_metal
#     ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
#     ni_metal = mol_ni - ni_sil

#     si_sil = (mol_o - ni_sil - fe_sil) / 2
#     si_metal = mol_si - si_sil

#     v_metal = bisection_search(g, root_bracket(g, mol_fe, mol_ni, ni_metal, mol_si, si_metal, mol_o, mol_v, v_metal, P_eq), mol_fe, 10e-7, mol_fe, mol_ni, mol_si, mol_o, P_eq)
#     v_sil = mol_v - v_metal
#     X_Si.append(si_metal)
#     X_FeO.append(fe_sil / (fe_sil + ni_sil + si_sil + v_sil))
#     X_Fe.append(fe_metal / (fe_metal + ni_metal + si_metal + v_metal))
#     X_V2O3.append(v_sil / (v_sil + ni_sil + fe_sil + si_sil))

# r_impactor = [i / 1e3 for i in r_impactor]


# fO2 = []
# for i in range(len(r_impactor)):
#     fO2.append(calculate_ln_o_iw_fugacity(X_FeO[i], X_Fe[i]))

# plt.plot(r_impactor, X_V2O3)
# plt.xlabel("Impactor radius (km)")
# plt.ylabel("X_V2O3")
# plt.title("X_V2O3 vs impactor radius (for a 1000 km radius planet)")
# plt.show()
