from fe_equilibrium import *

# knowns
r_planet = 100 * 1e3
r_impactor = 10 * 1e3
melt_factor = 10
k = 258838500680.3066 # using Mars' mass-g relationship

def calculate_g(M_p):
    """Calculate gravitational acceleration using g ~ M_p^(0.503) and M_Mars = 6.6e23 kg, g_Mars = 3.7 m / s"""
    return M_p**0.503 / k

def calculate_h(melt_vol, planet_size):
    # considering even melting across surface?
    """Returns the depth of mantle that is melted by an impactor."""
    return planet_size - (planet_size**3 - 3 / 4 / math.pi * melt_vol)**(1 / 3)

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
    return 4 / 3 * math.pi * (r_planet ** 3 - (r_planet - depth)**3) * rho_mantle

def convert_moles_to_mass(mol_element, element_molar_mass):
    """Returns the mass in kg of an element given its molar mass and number of moles."""
    return mol_element * element_molar_mass / 1000

def calculate_mantle_moles(element, mantle_depth, r_body):
    """Returns the amount of Si / V present in the mantle up to a specified depth."""
    if element == si or element == si_s:
        return convert_mass_to_moles(calculate_element_mass(calculate_shell_vol(mantle_depth, r_body) * rho_mantle, element), molar_mass_si)
    elif element == v or element == v_s:
        return convert_mass_to_moles(calculate_element_mass(calculate_shell_vol(mantle_depth, r_body) * rho_mantle, element), molar_mass_v)
    elif element == fe_s:
        return convert_mass_to_moles(calculate_element_mass(calculate_shell_vol(mantle_depth, r_body) * rho_mantle, element), molar_mass_fe)
    elif element == ni_s:
        return convert_mass_to_moles(calculate_element_mass(calculate_shell_vol(mantle_depth, r_body) * rho_mantle, element), molar_mass_ni)
    
fO2 = []
planet_size = []

# initial core radii
planet_core_radius = core_radius(r_planet)
c_impactor = core_radius(r_impactor)

# start with chondritic wt% - (Fe, Ni) in mantle
v_s = v
si_s = si
ni_s = 0
fe_s = 0
v_metal = 0

for i in range(10):
    planet_size.append(r_planet)
    h = calculate_h(melt_factor * calculate_vol(r_impactor), r_planet)
    g_acc = calculate_g(body_mass(r_planet, planet_core_radius))

    # calculate compounds already present in mantle
    mol_si = calculate_mantle_moles(si_s, h, r_planet)
    mol_v = calculate_mantle_moles(v_s, h, r_planet)
    mol_fe_s = calculate_mantle_moles(fe_s, h, r_planet)
    mol_ni_s = calculate_mantle_moles(ni_s, h, r_planet)

    # add compounds delivered. This gives total moles of each element in the mantle.
    mol_fe = calculate_total_element_moles(fe, molar_mass_fe, h, r_planet, r_impactor, c_impactor) + mol_fe_s
    mol_ni = calculate_total_element_moles(ni, molar_mass_ni, h, r_planet, r_impactor, c_impactor) + mol_ni_s
    mol_si += calculate_mantle_moles(si, r_impactor - c_impactor, r_impactor)
    mol_v += calculate_mantle_moles(v, r_impactor - c_impactor, r_impactor)
    mol_o = 2 * mol_si + mol_fe_s + mol_ni_s
    
    # equilibrium and mass balance
    P_eq = Peq(g_acc, h)
    fe_metal = bisection_search(f, root_bracket(f, mol_fe, mol_ni, mol_si, mol_o, mol_v, P_eq, v_metal), mol_fe, 10e-7, mol_fe, mol_ni, mol_si, mol_o, mol_v, P_eq, v_metal)
    fe_sil = mol_fe - fe_metal
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil) / 2
    si_metal = mol_si - si_sil
    v_metal = bisection_search(g, 0, root_bracket(g, mol_fe, mol_ni, mol_si, mol_o, mol_v, P_eq, fe_metal), 10e-7, mol_fe, mol_ni, mol_si, mol_o, mol_v, P_eq, fe_metal)
    v_sil = mol_v - v_metal

    # update mantle composition
    v_s = convert_moles_to_mass(v_sil, molar_mass_v) / body_mass(r_planet, planet_core_radius) * 100
    si_s = convert_moles_to_mass(si_sil, molar_mass_si) / body_mass(r_planet, planet_core_radius) * 100
    fe_s = convert_moles_to_mass(fe_sil, molar_mass_fe) / body_mass(r_planet, planet_core_radius) * 100
    print(fe_s)
    ni_s = convert_moles_to_mass(ni_sil, molar_mass_ni) / body_mass(r_planet, planet_core_radius) * 100

    # calculate fugacity
    X_FeO = fe_sil / (v_sil + fe_sil + ni_sil + si_sil)
    X_Fe = fe_metal / (fe_metal + v_metal + ni_metal + si_metal)
    fO2.append(calculate_ln_o_iw_fugacity(X_FeO, X_Fe))

    # increase planet size
    added_core_mass = convert_moles_to_mass(fe_metal, molar_mass_fe) + convert_moles_to_mass(ni_metal, molar_mass_ni) + convert_moles_to_mass(si_metal, molar_mass_si) + convert_moles_to_mass(v_metal, molar_mass_v)
    planet_core_radius += added_depth(planet_core_radius, added_core_mass, rho_core)
    r_planet += added_depth(r_planet, body_mass(r_impactor, c_impactor) - added_core_mass, rho_mantle)

plt.plot([r / 1e3 for r in planet_size], fO2)
plt.xlabel("Planet radius (km)")
plt.ylabel("ln(fO2) - IW")
plt.title("Fugacity vs. planet radius")
plt.show()
