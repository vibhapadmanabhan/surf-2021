from fe_equilibrium import *

# knowns
r_planet = 100 * 1e3
r_impactor = 10 * 1e3
melt_factor = 20
k = 258838500680.3066 # using Mars' mass-g relationship

def calculate_g(M_p):
    """Calculate gravitational acceleration using g ~ M_p^(0.503) and M_Mars = 6.6e23 kg, g_Mars = 3.7 m / s"""
    return M_p**0.503 / k

def calculate_h(melt_vol, planet_size):
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
    if element == si or element == si_s or element == si_s_impactor:
        return convert_mass_to_moles(calculate_element_mass(calculate_shell_vol(mantle_depth, r_body) * rho_mantle, element), molar_mass_si + molar_mass_o * 2)
    elif element == v or element == v_s or element == v_s_impactor:
        return convert_mass_to_moles(calculate_element_mass(calculate_shell_vol(mantle_depth, r_body) * rho_mantle, element), molar_mass_v * 2 + molar_mass_o * 3)
    elif element == fe_s:
        return convert_mass_to_moles(calculate_element_mass(calculate_shell_vol(mantle_depth, r_body) * rho_mantle, element), molar_mass_fe + molar_mass_o)
    elif element == ni_s:
        return convert_mass_to_moles(calculate_element_mass(calculate_shell_vol(mantle_depth, r_body) * rho_mantle, element), molar_mass_ni + molar_mass_o)
    elif element == mg_s or element == mg_s_impactor:
        return convert_mass_to_moles(calculate_element_mass(calculate_shell_vol(mantle_depth, r_body) * rho_mantle, element), molar_mass_mg + molar_mass_o)

def Teq(Peq):
    return 0.4 * 1661.2 * (Peq / 1.336 / 10**9 + 1)**(1 / 7.437)
    
# def scaled_wt_percent(metal, molar_mass_metal, half_ox_no, wt_percent, total_wt_percent):
#     return 

fO2 = []
planet_size = []

# initial core radii
planet_core_radius = core_radius(r_planet)
c_impactor = core_radius(r_impactor)

# start with chondritic wt% - (Fe, Ni) in mantle
fe_s = 8
mg = 9.54

# convert to weight percents of oxides and scale to make up entire mass
MgO = mg / molar_mass_mg * (molar_mass_mg + molar_mass_o)
SiO2 = si / molar_mass_si * (molar_mass_mg + 2 * molar_mass_o)
V2O3 = v / molar_mass_v * (molar_mass_v * 2 + molar_mass_o * 3)
FeO = fe_s / molar_mass_fe * (molar_mass_fe + molar_mass_o)
v_s = V2O3 * 100 / (V2O3 + SiO2 + MgO + FeO)
si_s = SiO2 * 100 / (V2O3 + SiO2 + FeO + MgO)
mg_s = MgO * 100 / (V2O3 + SiO2 + FeO + MgO)
fe_s = FeO * 100 / (V2O3 + SiO2 + FeO + MgO)
ni_s = 0
print(si_s, mg_s, v_s, fe_s)
# copy variables for these wt% being used in impactors
si_s_impactor = si_s
v_s_impactor = v_s
mg_s_impactor = mg_s
fe_s_impactor = fe_s

v_metal = 0
X_FeO = []
X_Si = []
X_Va = []

for i in range(1000):
    planet_size.append(r_planet)
    h = calculate_h(melt_factor * calculate_vol(r_impactor), r_planet)
    g_acc = calculate_g(body_mass(r_planet, planet_core_radius))
    # calculate compounds already present in mantle up till melted depth
    mol_si = calculate_mantle_moles(si_s, h, r_planet)
    mol_v = calculate_mantle_moles(v_s, h, r_planet)
    mol_fe_s = calculate_mantle_moles(fe_s, h, r_planet)
    mol_ni_s = calculate_mantle_moles(ni_s, h, r_planet)
    mol_mg = calculate_mantle_moles(mg_s, h, r_planet)

    # add compounds delivered. This gives total moles of each element in the mantle.
    # from impactor core
    mol_fe = calculate_total_element_moles(fe, molar_mass_fe, h, r_planet, r_impactor, c_impactor) + mol_fe_s
    mol_ni = calculate_total_element_moles(ni, molar_mass_ni, h, r_planet, r_impactor, c_impactor) + mol_ni_s
    # from impactor mantle
    mol_si += calculate_mantle_moles(si_s_impactor, r_impactor - c_impactor, r_impactor)
    mol_v += calculate_mantle_moles(v_s_impactor, r_impactor - c_impactor, r_impactor)
    mol_mg += calculate_mantle_moles(mg_s_impactor, r_impactor - c_impactor, r_impactor)
    mol_o = 2 * mol_si + mol_fe_s + mol_ni_s + mol_mg
    
    # increase h (account for volume of impactor's mantle)
    added_mantle_depth = added_depth(r_planet, body_mass(r_impactor, c_impactor), rho_mantle)
    h += added_mantle_depth

    # equilibrium and mass balance
    P_eq = Peq(g_acc, h) # * 10
    T_cmb = Teq(P_eq)
    fe_metal = bisection_search(f, root_bracket(f, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, v_metal), mol_fe, 10e-7, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, v_metal)
    fe_sil = mol_fe - fe_metal
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil - mol_mg) / 2
    si_metal = mol_si - si_sil
    # print(root_bracket(g, mol_fe, mol_ni, mol_si, mol_o, mol_v, P_eq, fe_metal))
    v_metal = bisection_search(g, 0, root_bracket(g, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, fe_metal), 10e-7, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, fe_metal)
    v_sil = mol_v - v_metal

    # calculate fugacity
    X_FeO = fe_sil / (v_sil + fe_sil + ni_sil + si_sil + mol_mg)
    X_Fe = fe_metal / (fe_metal + v_metal + ni_metal + si_metal)
    X_Si.append(si_metal / (fe_metal + v_metal + ni_metal + si_metal))
    fO2.append(calculate_ln_o_iw_fugacity(X_FeO, X_Fe))
    # print(X_FeO / (X_Fe + X_FeO))
    # increase planet size
    added_core_mass = convert_moles_to_mass(fe_metal, molar_mass_fe) + convert_moles_to_mass(ni_metal, molar_mass_ni) + convert_moles_to_mass(si_metal, molar_mass_si) + convert_moles_to_mass(v_metal, molar_mass_v)
    added_core_depth = added_depth(planet_core_radius, added_core_mass, rho_core)
    planet_core_radius += added_core_depth
    r_planet += added_mantle_depth + added_core_depth

    # update composition
    v_mantle_mass = convert_moles_to_mass(v_sil, molar_mass_v * 2 + molar_mass_o * 3)
    si_mantle_mass = convert_moles_to_mass(si_sil, molar_mass_si + molar_mass_o * 2)
    ni_mantle_mass = convert_moles_to_mass(ni_sil, molar_mass_ni + molar_mass_o)
    fe_mantle_mass = convert_moles_to_mass(fe_sil, molar_mass_fe + molar_mass_o)
    mg_mantle_mass = convert_moles_to_mass(mol_mg, molar_mass_mg + molar_mass_o)
    mantle_mass = v_mantle_mass + si_mantle_mass + ni_mantle_mass + fe_mantle_mass + mg_mantle_mass
    scale_factor = 1
    v_s = scale_factor * v_mantle_mass / mantle_mass * 100
    si_s = scale_factor * si_mantle_mass / mantle_mass * 100
    print(si_s)
    fe_s = scale_factor * fe_mantle_mass / mantle_mass * 100
    print(fe_s)
    ni_s = scale_factor * ni_mantle_mass / mantle_mass * 100
    mg_s = scale_factor * mg_mantle_mass / mantle_mass * 100
    print(mg_s + v_s + fe_s + ni_s + si_s)

plt.plot([r / 1e3 for r in planet_size], fO2)
plt.xlabel("Planet radius (km)")
plt.ylabel("ln(fO2) - IW")
plt.title("Fugacity vs. planet radius")
plt.show()
print(X_Si[-10:])
