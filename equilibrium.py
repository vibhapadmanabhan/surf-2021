import math
import numpy as np
import matplotlib.pyplot as plt
from atmosphere import Keq_FeO_H2O, Keq_FeO_CO2
from numba import jit, njit, float32

# constants
r_planet = 4500 * 1E3  # m
rho_mantle = 3000
rho_core = 8000

# P_eq = 0.7 * 25 # in GPa, assuming P_eq = 0.7 * P_cmb. CMB pressure depends on size of the planet.
T_eq = 4000
P_eq = 0.7 * 25  # GPa

# Mantle composition: chondritic, in wt%, ignoring trace elements and Ni, Fe.
si = 10.7
fe_s = 0
v = 0.00606
mg = 9.54

# Core composition, in wt%
fe = 85
ni = 15

# molar masses (g / mol)
molar_mass_fe = 55.8
molar_mass_ni = 58.7
molar_mass_si = 28
molar_mass_o = 16
molar_mass_v = 50.9415
molar_mass_mg = 24.3
molar_mass_h = 1
molar_mass_c = 12

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
    """Given an impactor and planet, calculates the total number of moles of Si, Ni, Fe, or FeO. Does not calculate moles of O."""
    if element == si or element == fe_s or element == v or element == mg:
        mass = (calculate_shell_vol(planet_mantle_depth, r_planet) +
                calculate_shell_vol(r_impactor - impactor_core_radius, r_impactor)) * rho_mantle
        element_mass = calculate_element_mass(mass, element)
        return convert_mass_to_moles(element_mass, element_molar_mass)
    mass = calculate_vol(impactor_core_radius) * rho_core
    element_mass = calculate_element_mass(mass, element)
    return convert_mass_to_moles(element_mass, element_molar_mass)


@jit(nopython=True)  # (float)
def calculate_kd(metal, T, P, a, b, c):
    """Calculates K_d using constants given in Rubie 2015 Table S1 and Rubie 2015 equation (3)."""
    if metal == "si":
        Kd = np.exp(a + b / T + (-155 * P + 2.26 * P**2 - 0.011 * P**3) / T)
    else:
        Kd = np.exp(a + b / T + c * P / T)
    return Kd


# actual distribution constants from Rubie
# actual_kd_si = calculate_kd("si", T_cmb, P_eq, 2.98, -15934, None)
# actual_kd_ni = calculate_kd("ni", T_cmb, P_eq, 1.06, 1553, -98)
# actual_kd_v = calculate_kd("v", T_cmb, P_eq, -0.48, -5063, 0)

@njit
def f(fe_metal, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, v_metal):
    """
    Returns the difference between the actual value of K_d(Si-Fe) and calculated value for K_d(Si-Fe). (Root of function will be found at the 
    correct value of K_d)
    """
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.06, 1553, -98)
    actual_kd_si = calculate_kd("si", T_eq, P_eq, 2.98, -15934, 0.)
    fe_sil = mol_fe - fe_metal
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil - mol_mg) / 2  # ignore mass of V oxides
    si_metal = mol_si - si_sil
    conc_si_metal = si_metal / (si_metal + fe_metal + ni_metal + v_metal)
    conc_si_sil = si_sil / (si_sil + fe_sil + ni_sil + (mol_v - v_metal) + mol_mg)
    conc_fe_metal = fe_metal / (fe_metal + ni_metal + si_metal + v_metal)
    conc_fe_sil = fe_sil / (fe_sil + ni_sil + si_sil + (mol_v - v_metal) + mol_mg)
    return conc_si_metal * conc_fe_sil**2 / conc_si_sil / conc_fe_metal**2 - actual_kd_si


@njit
def g(v_metal, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, fe_metal):
    """
    Returns the difference between the actual value of K_d(V-Fe) and calculated value for K_d(V-Fe). (Root of function will be found at the 
    correct value of K_d)
    """
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.06, 1553, -98)
    actual_kd_v = calculate_kd("v", T_eq, P_eq, -0.48, -5063, 0)
    fe_sil = mol_fe - fe_metal
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil - mol_mg) / 2  # ignore mass of V oxides
    si_metal = mol_si - si_sil
    conc_v_metal = v_metal / (si_metal + fe_metal + ni_metal + v_metal)
    conc_v_sil = (mol_v - v_metal) / (si_sil + fe_sil + ni_sil + (mol_v - v_metal) * 2 + mol_mg)
    conc_fe_metal = fe_metal / (fe_metal + ni_metal + si_metal + v_metal)
    conc_fe_sil = fe_sil / (fe_sil + ni_sil + si_sil + mol_v - v_metal + mol_mg)
    return conc_v_metal * conc_fe_sil**(3 / 2) / conc_v_sil / conc_fe_metal**(3 / 2) - actual_kd_v


@njit
def solve_H2O_atmosphere(mol_fe_mo, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, mol_h):
    """
    Returns the difference between the actual value of K_d(FeO - H2O) and calculated value for K_d(FeO - H2O). (Root of function will be found at the 
    correct value of K_d)

    Reaction is FeO + H2 <--> Fe + H2O
    """
    actual_kd = Keq_FeO_H2O(T_eq)
    mol_H2O = mol_o - mol_fe_mo
    mol_H2 = (mol_h - 2 * mol_H2O) / 2
    # print("mol h2", mol_H2)
    conc_fe_mo = mol_fe_mo / (mol_ni + mol_fe_mo + mol_si + mol_v + mol_mg)
    conc_H2O = mol_H2O / (mol_H2O + mol_H2)
    conc_H2 = mol_H2 / (mol_H2O + mol_H2)
    # print("kd error", conc_H2O / conc_H2 / conc_fe_mo - actual_kd)
    return conc_H2O / conc_H2 / conc_fe_mo - actual_kd

@njit
def solve_CO2_atmosphere(mol_fe_mo, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, mol_c):
    """
    Returns the difference between the actual value of K_d(FeO - CO2) and calculated value for K_d(FeO - CO2). (Root of function will be found at the 
    correct value of K_d).

    Reaction is FeO + CO <--> Fe + CO2
    """
    actual_kd = Keq_FeO_CO2(T_eq)
    mol_carbon_compounds = mol_o - mol_fe_mo
    conc_fe_mo = mol_fe_mo / (mol_ni + mol_fe_mo + mol_si + mol_v + mol_mg)
    CO2_CO = actual_kd * conc_fe_mo
    mol_CO = mol_carbon_compounds * 1 / (1 + CO2_CO)
    mol_CO2 = mol_carbon_compounds - mol_CO
    conc_CO = mol_CO / (mol_CO + mol_CO2)
    conc_CO2 = mol_CO2 / (mol_CO + mol_CO2)
    # print("kd error CO2", conc_CO2 / conc_CO / conc_fe_mo - actual_kd)
    return conc_CO2 / conc_CO / conc_fe_mo - actual_kd


@njit
def root_bracket(metal, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, aux):
    """Finds the lower/upper bound of an interval that can be used for the bisection search."""
    if metal == "fe":
        val = mol_fe
        FB = f(mol_fe, mol_fe, mol_ni, mol_si,
               mol_o, mol_v, mol_mg, P_eq, T_eq, aux)
        while val > 0 and FB * f(val, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, aux) > 0:
            val /= 1.00001  # multiply by small enough number for root to be found
        return val

    elif metal == "v":
        val = 1e-6
        FA = g(0, mol_fe, mol_ni, mol_si, mol_o,
               mol_v, mol_mg, P_eq, T_eq, aux)
        while val <= mol_v and FA * g(val, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, aux) > 0:
            val *= 1.0001
        return val

    elif metal == "mol_c_mo":
        val = mol_fe
        FB = solve_H2O_atmosphere(val, mol_fe, mol_ni, mol_si, mol_o,
               mol_v, mol_mg, P_eq, T_eq, aux)
        while val > 0 and FB * solve_CO2_atmosphere(val, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, aux) > 0:
            val /= 1.0000001
        # print("val", val)
        return val



def bisection_search(metal, a, b, eps, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, aux):
    """Aux refers to either fe_metal, v_metal, or mol_h, depending on which function kd(Si-Fe) or kd(V-Fe) is being searched for the root."""
    if (metal == "fe"):
        func = f
    elif (metal == "v"):
        func = g
    elif (metal == "mol_fe_mo"):
        func = solve_H2O_atmosphere
    elif (metal == "mol_c_mo"):
        func = solve_CO2_atmosphere
    while True:
        FA = func(a, mol_fe, mol_ni, mol_si, mol_o,
                  mol_v, mol_mg, P_eq, T_eq, aux)
        elem_metal = (a + b) / 2
        FP = func(elem_metal, mol_fe, mol_ni, mol_si,
                  mol_o, mol_v, mol_mg, P_eq, T_eq, aux)
        if np.abs(FP) <= eps:  # close enough to true value
            break
        if FA * FP > 0:
            a = elem_metal
        else:
            b = elem_metal
    return elem_metal


def calculate_ln_o_iw_fugacity(X_FeO, X_Fe):
    return 2 * math.log(X_FeO / X_Fe)

def calculate_fugacity(X_prod, X_reag, Keq):
    return (X_prod / X_reag / Keq)**2
