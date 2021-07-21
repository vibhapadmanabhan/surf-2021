from fe_equilibrium import *

# knowns
r_planet = 100 * 1e3
r_impactor = 20 * 1e3
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

def planet_mass(r_planet, core_radius):
    return 4 / 3 * math.pi * (core_radius**3 * rho_core + (r_planet - core_radius)**3 * rho_mantle)

def ocean_mass(r_planet, depth):
    return 4 / 3 * math.pi * (r_planet ** 3 - (r_planet - depth)**3) * rho_mantle

def convert_moles_to_mass(mol_element, element_molar_mass):
    """Returns the mass in kg of an element given its molar mass and number of moles."""
    return mol_element * element_molar_mass / 1000

def calculate_mantle_moles(element, mantle_depth, r_body):
    """Returns the amount of Si / V present in the mantle up to a specified depth."""
    if element == si or element == si_planet:
        return convert_mass_to_moles(calculate_element_mass(calculate_shell_vol(mantle_depth, r_body) * rho_mantle, element), molar_mass_si)
    elif element == v or element == v_planet:
        return convert_mass_to_moles(calculate_element_mass(calculate_shell_vol(mantle_depth, r_body) * rho_mantle, element), molar_mass_v)
    

fO2 = []
planet_size = []

planet_core_radius = core_radius(r_planet)
v_metal = 0
c_impactor = core_radius(r_impactor)
v_planet = v
si_planet = si


for i in range(10):
    planet_size.append(r_planet)
    h = calculate_h(melt_factor * calculate_vol(r_impactor), r_planet)
    g_acc = calculate_g(planet_mass(r_planet, planet_core_radius))
    print(g_acc)
    # mol_fe_s = calculate_total_element_moles(fe_s, molar_mass_fe + molar_mass_o, h, r_planet, r_impactor, c_impactor)

    # the planet melts. calculate moles of freed Si, V in planet's mantle
    mol_si = calculate_mantle_moles(si_planet, calculate_h(melt_factor * calculate_vol(r_impactor), r_planet), r_planet)

    # change the wt% of Si
    si_planet = convert_moles_to_mass(mol_si, molar_mass_si) / ocean_mass(r_planet, h) * 100

    # same for V
    mol_v = calculate_mantle_moles(v_planet, calculate_h(melt_factor * calculate_vol(r_impactor), r_planet), r_planet)

    # change the wt% of V
    v_planet = convert_moles_to_mass(mol_v, molar_mass_v) / ocean_mass(r_planet, h) * 100

    # delivered moles
    # Fe and Ni from core of impactor
    mol_fe = calculate_total_element_moles(fe, molar_mass_fe, h, r_planet, r_impactor, c_impactor)
    mol_ni = calculate_total_element_moles(ni, molar_mass_ni, h, r_planet, r_impactor, c_impactor)

    # Si and V from mantle of impactor
    mol_si += calculate_mantle_moles(si, r_impactor - c_impactor, r_impactor)

    mol_v += calculate_mantle_moles(v, r_impactor - c_impactor, r_impactor)

    # O calculated from SiO2 (ignore V oxides)
    mol_o = 2 * mol_si
    
    P_eq = Peq(g_acc, h)
    # print(P_eq)
    fe_metal = bisection_search(f, root_bracket(f, mol_fe, mol_ni, mol_si, mol_o, mol_v, P_eq, v_metal), mol_fe, 10e-7, mol_fe, mol_ni, mol_si, mol_o, mol_v, P_eq, v_metal)

    fe_sil = mol_fe - fe_metal
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil) / 2
    si_metal = mol_si - si_sil

    v_metal = bisection_search(g, 0, root_bracket(g, mol_fe, mol_ni, mol_si, mol_o, mol_v, P_eq, fe_metal), 10e-7, mol_fe, mol_ni, mol_si, mol_o, mol_v, P_eq, fe_metal)

    v_sil = mol_v - v_metal

    X_FeO = fe_sil / (v_sil + fe_sil + ni_sil + si_sil)
    X_Fe = fe_metal / (fe_metal + v_metal + ni_metal + si_metal)

    fO2.append(calculate_ln_o_iw_fugacity(X_FeO, X_Fe))

    added_core_mass = convert_moles_to_mass(fe_metal, molar_mass_fe) + convert_moles_to_mass(ni_metal, molar_mass_ni) + convert_moles_to_mass(si_metal, molar_mass_si) + convert_moles_to_mass(v_metal, molar_mass_v)

    planet_core_radius += added_depth(planet_core_radius, added_core_mass, rho_core)

    r_planet += added_depth(r_planet, calculate_shell_vol(r_impactor - core_radius(r_impactor), r_impactor) * rho_mantle, rho_mantle)

plt.plot([r / 1e3 for r in planet_size], fO2)
plt.xlabel("Planet radius (km)")
plt.ylabel("ln(fO2) - IW")
plt.title("Fugacity vs. planet radius")
plt.show()
