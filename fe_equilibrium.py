import math

# constants 
r_planet = 1000 * 10**3
r_impactor = 100 * 10**3
rho_mantle = 3000 
rho_core = 8000
P_cmb = 136 * 0.7 # in GPa, assuming 70% of CMB pressure
T_cmb = 4000

# Mantle composition: chondritic, in wt%, ignoring trace elements.
si = 10.7

# Core composition, in wt%
fe = 85
ni = 15

planet_mantle_depth = 300 * 10**3  # given
impactor_core_radius = 53.1329  * 10**3  # assuming core/mantle mass ratio to be the same as present day Earth's
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
kd_si = calculate_kd("si", T_cmb, P_cmb, 2.98, -15934, None)
kd_ni = calculate_kd("ni", T_cmb, P_cmb, 1.06, 1553, -98)
print("Actual K_d for Ni-Fe: ", kd_ni)
print("Actual K_d for Si-Fe: ", kd_si)


# assume some initial values for mols of x in the metal phase (core)
init_mol_fe_metal = 10**17
init_mol_si_metal = 6 * 10**20

# mass balance
def calculate_vals(mol_fe_metal, mol_si_metal, mol_si, mol_ni, mol_o, mol_fe):
    mol_FeO = mol_fe - mol_fe_metal
    mol_SiO2 = mol_si - mol_si_metal
    mol_NiO = mol_o - 2 * mol_SiO2 - mol_FeO
    mol_ni_metal = mol_ni - mol_NiO
    return ((mol_fe_metal, mol_FeO, mol_ni_metal, mol_NiO, mol_si_metal, mol_SiO2))

mol_fe_metal, mol_FeO, mol_ni_metal, mol_NiO, mol_si_metal, mol_SiO2 = calculate_vals(init_mol_fe_metal, init_mol_si_metal, mol_si, mol_ni, mol_o, mol_fe)

# print(mol_fe_metal, mol_FeO, mol_ni_metal, mol_NiO, mol_si_metal, mol_SiO2)
conc_fe_metal = mol_fe_metal / (mol_fe_metal + mol_si_metal + mol_ni_metal)
conc_FeO = mol_FeO / (mol_FeO + mol_NiO + mol_SiO2)
conc_ni_metal = mol_ni_metal / (mol_fe_metal + mol_si_metal + mol_ni_metal)
conc_NiO = mol_NiO / (mol_FeO + mol_NiO + mol_SiO2)
conc_si_metal = mol_si_metal / (mol_fe_metal + mol_si_metal + mol_ni_metal)
conc_SiO2 = mol_SiO2 / (mol_FeO + mol_NiO + mol_SiO2)

print("Calculated value for Ni-Fe: ", conc_ni_metal * conc_FeO / conc_NiO / conc_fe_metal)
print("Calculated value for Si-Fe: ", conc_si_metal * conc_FeO / conc_SiO2 / conc_fe_metal)
