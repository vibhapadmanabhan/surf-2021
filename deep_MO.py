from equilibrium import *
from growth import *
"""This file models the deep magma ocean concept. The magma ocean is the entire mantle."""

planet_core_radius = core_radius(r_planet)

# in the planet
mol_fe = 1000 * fe_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_fe)
mol_ni = 1000 * ni_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_ni)
mol_si = 1000 * si_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_si)
mol_mg = 1000 * mg_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_mg)
mol_v = 1000 * v_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_v)
mol_o = 1000 * o_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_o)

mols_si_c = 1000 * si * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_si
mols_v_c = 1000 * v * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_v
mols_ni_c = 1000 * ni * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_ni
mols_fe_c = 1000 * fe * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_fe

v_metal = 0
# fe_s = [i for i in range(20)]
# si_s = [(54 - 0.00606 - i / 2) / (molar_mass_si + molar_mass_o * 2) * molar_mass_si for i in fe_s]
# mg_s = [(46 - i / 2) / (molar_mass_mg + molar_mass_o) * molar_mass_mg for i in fe_s]
# o_s = []
# for i in range(len(fe_s)):
#     o_s.append(100 - fe_s[i] - mg_s[i] - si_s[i] - v_s)

j = 0
fe_s = 17.48
si_s = (54 - 0.00606 - fe_s / 2) / (molar_mass_si + molar_mass_o * 2) * molar_mass_si
mg_s = (46 - fe_s / 2) / (molar_mass_mg + molar_mass_o) * molar_mass_mg
o_s = 100 - fe_s - mg_s - si_s - v_s

while r_planet <= 1000e3:
    print(j)
    planet_size.append(r_planet)
    X_FeO_impactor.append(fe_s)
    # growing in size
    r_impactor = r_planet / 5
    impactor_size.append(r_impactor)
    c_impactor = core_radius(r_impactor)
    delivered_fe = 1000 * fe * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_fe
    delivered_ni = 1000 * ni * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_ni
    delivered_fe_s = fe_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) * 1000 / molar_mass_fe
    delivered_si = 1000 * si_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_si
    delivered_mg = 1000 * mg_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_mg
    delivered_v = 1000 * v_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_v
    delivered_o = 1000 * o_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_o
    h = r_planet - planet_core_radius

    planet_mass = convert_moles_to_mass(mol_fe, molar_mass_fe) + convert_moles_to_mass(mol_ni, molar_mass_ni) + convert_moles_to_mass(mol_si, molar_mass_si) + convert_moles_to_mass(mol_v, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg) + convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_v_c, molar_mass_v)
    g_acc = calculate_g(planet_mass)
    gravity.append(g_acc)
    
    # add compounds delivered. This gives total moles of each element in the mantle.
    mol_ni += delivered_ni
    # from impactor mantle
    mol_si += delivered_si
    mol_v += delivered_v
    mol_mg += delivered_mg
    mol_fe += delivered_fe_s + delivered_fe
    mol_o += delivered_o

    # equilibrium and mass balance
    P_eq = Peq(g_acc, h) # assume CMB pressure
    pressure.append(P_eq)
    T_eq = Teq(P_eq * 1e9)
    temperature.append(T_eq)
    # small error due to exclusion of delta h from volume of impactor mantle (does not change results)
    fe_metal = bisection_search("fe", root_bracket("fe", mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, v_metal), mol_fe, 10e-12, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, v_metal)
    fe_sil = mol_fe - fe_metal
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.06, 1553, -98)
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil - mol_mg) / 2
    si_metal = mol_si - si_sil
    v_metal = bisection_search("v", 0, root_bracket('v', mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, fe_metal), 10e-12, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, fe_metal)
    v_sil = mol_v - v_metal
    # add moles in metal phase to total moles of each element in the core
    mols_fe_c += fe_metal
    mols_ni_c += ni_metal
    mols_si_c += si_metal
    mols_v_c += v_metal

    # calculate new mantle and core mass
    new_mantle_mass = convert_moles_to_mass(fe_sil, molar_mass_fe) + convert_moles_to_mass(ni_sil, molar_mass_ni) + convert_moles_to_mass(si_sil, molar_mass_si) + convert_moles_to_mass(v_sil, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg)
    new_core_mass = convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_v_c, molar_mass_v)

    # increase planet size
    planet_core_radius = sphere_radius(new_core_mass, rho_core)
    new_mantle_depth = shell_width(new_mantle_mass, rho_mantle, planet_core_radius)
    mantle_depth.append(new_mantle_depth)
    r_planet = planet_core_radius + new_mantle_depth
    # update total moles in mantle of elements involved in equilibrium based on number of moles left in the silicate phase
    mol_fe = fe_sil
    mol_ni = ni_sil
    mol_si = si_sil
    mol_v = v_sil

    # change impactor composition
    fe_s -= 0.02
    si_s = (54 - 0.00606 - fe_s / 2) / (molar_mass_si + molar_mass_o * 2) * molar_mass_si
    mg_s = (46 - fe_s / 2) / (molar_mass_mg + molar_mass_o) * molar_mass_mg
    o_s = 100 - fe_s - mg_s - si_s - v_s

    # calculate fugacity
    X_FeO.append(mol_fe / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    X_Fe.append(mols_fe_c / (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c))
    fO2.append(calculate_ln_o_iw_fugacity(X_FeO[-1], X_Fe[-1]))

    # track other values
    X_Si.append(mols_si_c / (mols_fe_c + mols_v_c + mols_ni_c + mols_si_c))
    X_Mg.append(mol_mg / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    X_NiO.append(mol_ni / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    X_Ni.append(mols_ni_c / (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c))
    X_SiO2.append(mol_si / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    X_Va.append(mols_v_c / (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c))
    X_VO.append(mol_v * 2 / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    j += 1

save_data(X_Fe, X_Si, X_Ni, X_Va, X_FeO, X_SiO2, X_NiO, X_Mg, X_VO, X_FeO_impactor, gravity, pressure, temperature, planet_size, impactor_size, mantle_depth, fO2, "./data/deep-MO/heterogeneous_accretion_decreasing_FeO.txt")
