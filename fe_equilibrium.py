import math
from scipy.optimize import fsolve

# constants
r_planet = 1000 * 10E3
r_impactor = 100 * 10E3
rho_mantle = 3000 
rho_core = 8000
P_cmb = 0.7 * 25 # in GPa, assuming 70% of CMB pressure. CMB pressure depends on size of the planet.
T_cmb = 4000

# Mantle composition: chondritic, in wt%, ignoring trace elements.
si = 10.7

# Core composition, in wt%
fe = 85
ni = 15

planet_mantle_depth = 300 * 10E3  # given
impactor_core_radius = 53.1329  * 10E3  # assuming core/mantle mass ratio to be the same as present day Earth's
impactor_mantle_depth = r_impactor - impactor_core_radius

def calculate_shell_vol(mantle_depth, r_obj):
    return 4 / 3 * math.pi * (r_obj**3 - (r_obj - mantle_depth)**3)

def calculate_vol(r_obj):
    return 4 / 3 * math.pi * r_obj**3

def calculate_element_mass(section_mass, element_wt_percent):
    return section_mass * element_wt_percent * 0.01

def convert_mass_to_moles(compound_mass, compound_molar_mass):
    return compound_mass / compound_molar_mass

impactor_mantle_mass = calculate_shell_vol(impactor_mantle_depth, r_impactor) * rho_mantle
planet_mantle_mass = calculate_shell_vol(planet_mantle_depth, r_planet) * rho_mantle
impactor_core_mass = calculate_vol(impactor_core_radius) * rho_core

fe_mass = calculate_element_mass(impactor_core_mass, fe)
ni_mass = calculate_element_mass(impactor_core_mass, ni)

# total amounts of elements
mol_fe = convert_mass_to_moles(fe_mass, 55.8)
mol_ni = convert_mass_to_moles(ni_mass, 58.7)
mol_si = convert_mass_to_moles(calculate_element_mass(planet_mantle_mass, si) + calculate_element_mass(impactor_mantle_mass, si), 28)
mol_o = 2 * mol_si  # from SiO2

def calculate_kd(metal, T, P, a, b, c):
    if metal == "si":
        return math.exp(a + b / T + (-155 * P + 2.26 * P**2 - 0.011 * P**3) / T)
    else:
        return 10**(a + b / T + c * P / T)

# actual distribution constants from Rubie
actual_kd_si = calculate_kd("si", T_cmb, P_cmb, 2.98, -15934, None)
actual_kd_ni = calculate_kd("ni", T_cmb, P_cmb, 1.06, 1553, -98)
print("Actual k_d", actual_kd_ni, actual_kd_si)
# mass balance
def calculate_vals(mol_fe_metal, mol_si_metal, mol_si, mol_ni, mol_o, mol_fe):
    """
    Returns the number of moles of (FeO, Ni, NiO, SiO2) in a tuple given initial values for Fe and Si.
    """
    mol_FeO = mol_fe - mol_fe_metal
    mol_SiO2 = mol_si - mol_si_metal
    mol_NiO = mol_o - 2 * mol_SiO2 - mol_FeO
    mol_ni_metal = mol_ni - mol_NiO
    return ((mol_FeO, mol_ni_metal, mol_NiO, mol_SiO2))

def kd_ni(mol_fe_metal, mol_FeO, mol_ni_metal, mol_NiO, mol_si_metal, mol_SiO2):
    """
    Calculates K_d for the reaction NiO + Fe <--> FeO + Ni.
    """
    conc_fe_metal = mol_fe_metal / (mol_fe_metal + mol_si_metal + mol_ni_metal)
    conc_FeO = mol_FeO / (mol_FeO + mol_NiO + mol_SiO2)
    conc_ni_metal = mol_ni_metal / (mol_fe_metal + mol_si_metal + mol_ni_metal)
    conc_NiO = mol_NiO / (mol_FeO + mol_NiO + mol_SiO2)
    return conc_ni_metal * conc_FeO / conc_NiO / conc_fe_metal

def kd_si(mol_fe_metal, mol_FeO, mol_ni_metal, mol_NiO, mol_si_metal, mol_SiO2):
    """
    Calculates K_d for the reaction # SiO2 + 2Fe <--> 2FeO + Si.
    """
    conc_fe_metal = mol_fe_metal / (mol_fe_metal + mol_si_metal + mol_ni_metal)
    conc_FeO = mol_FeO / (mol_FeO + mol_NiO + mol_SiO2)
    conc_si_metal = mol_si_metal / (mol_fe_metal + mol_si_metal + mol_ni_metal)
    conc_SiO2 = mol_SiO2 / (mol_FeO + mol_NiO + mol_SiO2)
    return conc_si_metal * conc_FeO**2 / conc_SiO2 / conc_fe_metal**2

# initial guesses for Fe and Si
mol_fe_metal = mol_fe - 10E17
mol_si_metal = 10E15

# calculate other initial values
mol_FeO, mol_ni_metal, mol_NiO, mol_SiO2 = calculate_vals(mol_fe, 0, mol_si, mol_ni, mol_o, mol_fe)

# print initial values
print("Initial values:")
print(mol_fe_metal, mol_FeO, mol_ni_metal, mol_NiO, mol_si_metal, mol_SiO2)


kd_ni_calc = kd_ni(mol_fe_metal, mol_FeO, mol_ni_metal, mol_NiO, mol_si_metal, mol_SiO2)

# use Le Chatelier's principle to converge on correct values
while kd_ni_calc < actual_kd_ni:
    # print("loop")
    mol_fe_metal += 10E17
    mol_FeO, mol_ni_metal, mol_NiO, mol_SiO2 = calculate_vals(mol_fe_metal, mol_si_metal, mol_si, mol_ni, mol_o, mol_fe)
    kd_si_calc = kd_si(mol_fe_metal, mol_FeO, mol_ni_metal, mol_NiO, mol_si_metal, mol_SiO2)
    if kd_si_calc < actual_kd_si:
        mol_si_metal += 10E17
        mol_FeO, mol_ni_metal, mol_NiO, mol_SiO2 = calculate_vals(mol_fe_metal, mol_si_metal, mol_si, mol_ni, mol_o, mol_fe)
    kd_ni_calc = kd_ni(mol_fe_metal, mol_FeO, mol_ni_metal, mol_NiO, mol_si_metal, mol_SiO2)


print("final values: ", mol_fe_metal, mol_si_metal)
print("calculated K_d: ", kd_ni(mol_fe_metal, mol_FeO, mol_ni_metal, mol_NiO, mol_si_metal, mol_SiO2), kd_si(mol_fe_metal, mol_FeO, mol_ni_metal, mol_NiO, mol_si_metal, mol_SiO2))
# print("Calculated value for Ni-Fe: ", kd_ni(mol_fe_metal, mol_FeO, mol_ni_metal, mol_NiO, mol_si_metal, mol_SiO2))

# print("Calculated value for Si-Fe: ", kd_si(mol_fe_metal, mol_FeO, mol_ni_metal, mol_NiO, mol_si_metal, mol_SiO2))


# x[0] is Fe_m, x[1] is FeO, x[2] is Si_m, x[3] is SiO, x[4] is Ni_m, x[5] is NiO
# def func(x, mol_fe, mol_ni, mol_si, mol_o, kd_ni, kd_si):
#    return [x[0] + x[1] - mol_fe, x[2] + x[3] - mol_si, x[4] + x[5] - mol_ni, x[1] + x[3] + x[5] - mol_o, 
#    x[4] / (x[4] + x[0] + x[2]) * (x[1] / (x[1] + x[3] + x[5])) / (x[5] / (x[1] + x[2] + x[3])) / (x[0] / (x[2] + x[4] + x[0])) - kd_ni, 
#    x[2] / (x[4] + x[0] + x[2]) * (x[1] / (x[1] + x[3] + x[5]))**2 / (x[3] / (x[1] + x[2] + x[3])) / (x[0] / (x[2] + x[4] + x[0]))**2 -  #kd_si]

#root = fsolve(func, [10E19, 10E18, 10E17, 10E21, 10E20, 10E10], (mol_fe, mol_ni, mol_si, mol_o, kd_ni, kd_si))
#print(root)
