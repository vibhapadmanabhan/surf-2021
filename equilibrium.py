import math
import numpy as np
import matplotlib.pyplot as plt
from atmosphere import *
from calculate_initial_parameters import *
from volumes import *
from masses import *
from densities import *
from planetary_embryo_initial_parameters import *
from conversions import *
from initial_composition import *

planet_core_radius = core_radius(r_planet)
mantle_mass = calculate_shell_vol(r_planet - planet_core_radius, r_planet) * rho_mantle
core_mass = calculate_vol(planet_core_radius) * rho_core
P_eq = Peq(calculate_g(mantle_mass + core_mass), planet_mantle_depth) # GPa
T_eq = Teq(P_eq * 1e9)


def calculate_impactor_core_radius(impactor_radius):
    """
    Returns the radius of the core of a body in meters using Earth's core-mantle mass ratio.
    """
    return (0.15 * impactor_radius**3)**(1 / 3)


def calculate_element_mass(section_mass, element_wt_percent):
    """Returns the mass of an element in g in the mantle/core given the wt% of the element. """
    return section_mass * element_wt_percent * 0.01 * 1000


def calculate_total_element_moles(element, element_molar_mass, planet_mantle_depth, r_planet, r_impactor, impactor_core_radius):
    """Given an impactor and planet, calculates the total number of moles of Si, Ni, Fe, or FeO. Does not calculate moles of O."""
    if element == si_s or element == fe_s or element == v_s or element == mg_s or element == o_s:
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
    correct value of K_d). The equilibrium reactions occurring are:

    2 Fe + SiO2 <--> 2 FeO + Si
    2 Ni + SiO2 <--> 2 NiO + Si
    """
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.06, 1553, -98)
    actual_kd_si = calculate_kd("si", T_eq, P_eq, 2.98, -15934, 0.)
    fe_sil = mol_fe - fe_metal
    # print("total mol fe, total mol fe metal", mol_fe, fe_metal)
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil - mol_mg) / 2  # ignore mass of V oxides
    si_metal = mol_si - si_sil
    conc_si_metal = si_metal / (si_metal + fe_metal + ni_metal + v_metal)
    conc_si_sil = si_sil / (si_sil + fe_sil + ni_sil + (mol_v - v_metal) + mol_mg)
    conc_fe_metal = fe_metal / (fe_metal + ni_metal + si_metal + v_metal)
    conc_fe_sil = fe_sil / (fe_sil + ni_sil + si_sil + (mol_v - v_metal) + mol_mg)
    # print("conc si metal    conc fe sil     conc si sil     conc fe metal", conc_si_metal, conc_fe_sil, conc_si_sil, conc_fe_metal)
    return conc_si_metal * conc_fe_sil**2 / conc_si_sil / conc_fe_metal**2 - actual_kd_si

# max value for fe metal = 
# fe_sil = mol_o - ni_sil - mol_mg - mol_si * 2

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
    # print("kd difference v  ", conc_v_metal * conc_fe_sil**(3 / 2) / conc_v_sil / conc_fe_metal**(3 / 2) - actual_kd_v)
    return conc_v_metal * conc_fe_sil**(3 / 2) / conc_v_sil / conc_fe_metal**(3 / 2) - actual_kd_v
    

def solve_initial_atmosphere(mol_H2O, mol_o_atmos, mol_h, mol_c, fO2, T_eq):
    """
    Solves for initial amounts of all species, given fO2, where equilibria are governed by the following reactions:

    CO + 0.5 O2 <--> CO2
    CH4 + 2 O2 <--> CO2 + 2H2O
    H2 + 0.5 O2 <--> H2O
    """
    mol_H2 = mol_H2O / H2O_H2ratio(fO2, T_eq)
    mol_CH4 = (mol_h - 2 * mol_H2 - 2 * mol_H2O) / 4
    mol_CO2 = mol_o_atmos - mol_c - mol_H2O + mol_CH4
    mol_CO = mol_CO2 / CO2_COratio(fO2, T_eq)
    mol_volatiles = mol_H2O + mol_H2 + mol_CO + mol_CO2 + mol_CH4
    conc_CO2 = mol_CO2 / mol_volatiles
    conc_H2O = mol_H2O / mol_volatiles
    conc_CH4 = mol_CH4 / mol_volatiles
    return conc_CO2 * conc_H2O**2 / conc_CH4 - CO2H2O_CH4ratio(fO2, T_eq)


@njit
def solve_atmosphere(mol_fe_mo, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, mol_c, T_eq, mol_h):
    """
    By finding the difference between the actual Kd(FeO - CH4) and calculated Kd(FeO - CH4), solves for FeO, Fe, CO, CO2, H2, H2O molar amounts during re-equilibrium governed by the following equations:

    FeO + CO <--> Fe + CO2
    FeO + CH4 <--> Fe + CO + 2 H2
    FeO + H2 <--> Fe + H2O

    CH4 is ignored due to its low concentration.
    """
    kd_CO2 = Keq_FeO_CO2(T_eq)
    # kd_CH4 = Keq_FeO_CH4(T_eq)
    kd_H2O = Keq_FeO_H2O(T_eq)
    mol_fe_metal = mol_fe - mol_fe_mo
    conc_fe_mo = mol_fe_mo / (mol_ni + mol_fe_mo + mol_si + mol_mg + mol_v)
    H2O_to_H2 = kd_H2O * conc_fe_mo
    # ignore CH4 here because its concentration is low
    mol_H2 = 1 / (H2O_to_H2 + 1) * mol_h / 2
    mol_H2O = (mol_h - 2 * mol_H2) / 2
    CO2_to_CO = kd_CO2 * conc_fe_mo
    mol_CO = 1 / (CO2_to_CO + 1) * mol_c
    mol_CO2 = mol_c - mol_CO
    # mol_CH4 = mol_c - mol_CO - mol_CO2
    mol_volatiles = mol_H2O + mol_H2 + mol_CO + mol_CO2
    # conc_CH4 = mol_CH4 / mol_volatiles
    conc_CO = mol_CO / mol_volatiles
    conc_CO2 = mol_CO2 / mol_volatiles
    #print("co, co2, h2, h2o, mol fe", mol_CO, mol_CO2, mol_H2, mol_H2O, mol_fe)
    #print("error", (mol_o - mol_CO - mol_CO2 * 2 - mol_H2O - mol_fe_mo) / mol_o)
    return (mol_o - mol_CO - mol_CO2 * 2 - mol_H2O - mol_fe_mo) / mol_o

@njit
def root_bracket(metal, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, aux):
    """Finds the lower/upper bound of an interval that can be used for the bisection search."""
    if metal == "fe":
        val = 1
        FB = f(mol_fe, mol_fe, mol_ni, mol_si,
               mol_o, mol_v, mol_mg, P_eq, T_eq, aux)
        while val < mol_fe and FB * f(val, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, aux) > 0:
            val *= 1.00000001  # multiply by small enough number for root to be found
        # print("lower bound", val)
        return val

    elif metal == "v":
        val = 1e-6
        FA = g(0, mol_fe, mol_ni, mol_si, mol_o,
               mol_v, mol_mg, P_eq, T_eq, aux)
        while val < mol_v and FA * g(val, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, aux) > 0:
            val *= 1.001
        return val


def bisection_search(metal, a, b, eps, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, aux):
    """Aux refers to either fe_metal, v_metal, or mol_h, depending on which function kd(Si-Fe) or kd(V-Fe) is being searched for the root."""
    if (metal == "fe"):
        func = f
    elif (metal == "v"):
        func = g
    elif (metal == "mol_fe_mo"):
        func = solve_atmosphere
    elif (metal == "atmos_initial"):
        func = solve_initial_atmosphere

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


def bisection_search_atmosphere(a, b, eps, mol_o_atmos, mol_h, mol_c, fO2, T_eq):
    """Solves for the molar amount of H2O in the atmosphere after equilibrium, given the fO2."""
    while True:
        FA = solve_initial_atmosphere(a, mol_o_atmos, mol_h, mol_c, fO2, T_eq)
        elem_metal = (a + b) / 2
        FP = solve_initial_atmosphere(elem_metal, mol_o_atmos, mol_h, mol_c, fO2, T_eq)

        if np.abs(FP) <= eps:  # close enough to true value
            break
        if FA * FP > 0:
            a = elem_metal
        else:
            b = elem_metal
    return elem_metal

def calculate_ln_o_iw_fugacity(X_FeO, X_Fe):
    """Calculates fO2 - IW"""
    return 2 * math.log(X_FeO / X_Fe)

def calculate_fugacity(X_prod, X_reag, Keq):
    """Calculates fO2"""
    return (X_prod / X_reag / Keq)**2

def calc_metal(mol_fe, mol_ni, mol_si, mol_mg, T_eq, P_eq):
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.6, 1553, -98)
    actual_kd_si = calculate_kd("si", T_eq, P_eq, 2.98, -15934, 0)

    mol_oxide = mol_fe + mol_ni + 2 * mol_si + mol_mg

    r_NiO = mol_ni / mol_oxide
    r_FeO = mol_fe / mol_oxide
    r_SiO2 = mol_si / mol_oxide
    r_NiFe = actual_kd_ni * r_NiO / r_FeO
    r_SiFe2 = actual_kd_si * r_SiO2 / (r_FeO * r_FeO)

