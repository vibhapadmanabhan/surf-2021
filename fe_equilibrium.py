import math

# constants 
r_planet = 1000 * 10**3
r_impactor = 100 * 10**3
rho_mantle = 3000 
rho_core = 8000

# Mantle composition: all chondritic, in wt%
si = 10.7
mg = 9.54
al = 1.02
ca = 1.11
co = 0.0513
v = 0.0061
cr = 0.2623

# Core composition, in wt%
fe = 85
ni = 15

planet_mantle_depth = 300 * 10**3  # given
impactor_core_radius = 53.1329  * 10**3  ## assuming core/mantle mass ratio to be the same as present day Earth's
impactor_mantle_depth = r_impactor - impactor_core_radius

def calculate_shell_vol(mantle_depth, r_obj):
    return 4 / 3 * math.pi * (r_obj**3 - (r_obj - mantle_depth)**3)

def calculate_vol(r_obj):
    return 4 / 3 * math.pi * r_obj**3

def calculate_element_mass(section_mass, element_wt_percent):
    return section_mass * element_wt_percent * 0.01

impactor_mantle_mass = calculate_shell_vol(impactor_mantle_depth, r_impactor) * rho_mantle
planet_mantle_mass = calculate_shell_vol(planet_mantle_depth, r_planet) * rho_mantle

impactor_core_mass = calculate_vol(impactor_core_radius) * rho_core
fe_mass = calculate_element_mass(impactor_core_mass, fe)
# ni_mass = calculate_element_mass(impactor_core_mass, ni)

si_mass = calculate_element_mass(planet_mantle_mass, si) + calculate_element_mass(impactor_mantle_mass, si)
mg_mass = calculate_element_mass(planet_mantle_mass, si)
al_mass = calculate_element_mass(planet_mantle_mass, si)
co_mass = calculate_element_mass(planet_mantle_mass, si)
v_mass = calculate_element_mass(planet_mantle_mass, si)
cr_mass = calculate_element_mass(planet_mantle_mass, si)

def calculate_kd(metal, T, P, a, b, c):
    if metal == "si":
        return math.exp(a + b / T + (-155 * P + 2.26 * P**2 - 0.011 * P**3) / T)
    else:
        return 10**(a + b / T + c * P / T)

