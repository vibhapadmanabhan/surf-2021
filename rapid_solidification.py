from equilibrium import *
from growth import *
"""This file models the rapid solidification theory. The magma ocean solidifies completely before the next impactor arrives."""

planet_core_radius = core_radius(r_planet)
melt_factor = 10
h = r_planet - planet_core_radius
# in the planet
mol_fe = 1000 * fe_s * 0.01 * ocean_mass(r_planet, h) / molar_mass_fe
mol_ni = 1000 * ni_s * 0.01 * ocean_mass(r_planet, h) / molar_mass_ni
mol_si = 1000 * si_s * 0.01 * ocean_mass(r_planet, h) / molar_mass_si
mol_mg = 1000 * mg_s * 0.01 * ocean_mass(r_planet, h) / molar_mass_mg
mol_v = 1000 * v_s * 0.01 * ocean_mass(r_planet, h) / molar_mass_v
mol_o = 1000 * o_s * 0.01 * ocean_mass(r_planet, h) / molar_mass_o

mols_si_c = 1000 * si * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_si
mols_v_c = 1000 * v * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_v
mols_ni_c = 1000 * ni * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_ni
mols_fe_c = 1000 * fe * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_fe

v_metal = 0
j = 0
fe_s = 17.52
si_s = (54 - 0.00606 - fe_s / 2) / (molar_mass_si + molar_mass_o * 2) * molar_mass_si
mg_s = (46 - fe_s / 2) / (molar_mass_mg + molar_mass_o) * molar_mass_mg
o_s = 100 - fe_s - mg_s - si_s - v_s

while r_planet <= 1000e3:
    print(j)
    planet_size.append(r_planet)
    r_impactor = r_planet / 5
    impactor_size.append(r_impactor)
    c_impactor = core_radius(r_impactor)
    X_FeO_impactor.append(fe_s)

    delivered_fe = 1000 * fe * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_fe
    delivered_ni = 1000 * ni * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_ni

    delivered_fe_s = fe_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) * 1000 / molar_mass_fe
    delivered_si = 1000 * si_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_si
    delivered_mg = 1000 * mg_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_mg
    delivered_v = 1000 * v_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_v
    delivered_o = 1000 * o_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_o
    h = calculate_h(melt_factor * calculate_vol(r_impactor), r_planet)
    # h_frac is the fraction of the mantle that is molten (magma ocean)
    h_frac = calculate_shell_vol(h, r_planet) / calculate_shell_vol(r_planet - planet_core_radius, r_planet)
    planet_mass = convert_moles_to_mass(mol_fe, molar_mass_fe) + convert_moles_to_mass(mol_ni, molar_mass_ni) + convert_moles_to_mass(mol_si, molar_mass_si) + convert_moles_to_mass(mol_v, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg) + convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_v_c, molar_mass_v)
    g_acc = calculate_g(planet_mass)
    gravity.append(g_acc)

    impactor_mass = convert_moles_to_mass(delivered_fe, molar_mass_fe) + convert_moles_to_mass(delivered_ni, molar_mass_ni) + convert_moles_to_mass(delivered_si, molar_mass_si) + convert_moles_to_mass(delivered_v, molar_mass_v) + convert_moles_to_mass(delivered_o, molar_mass_o) + convert_moles_to_mass(delivered_mg, molar_mass_mg) + convert_moles_to_mass(delivered_fe_s, molar_mass_fe) 

    # calculate compounds already present in mantle up till melted depth assuming homogeneous mantle
    mol_si_MO = h_frac * mol_si
    mol_v_MO = h_frac * mol_v
    mol_FeO_MO = h_frac * mol_fe
    mol_ni_MO = h_frac * mol_ni
    mol_mg_MO = h_frac * mol_mg
    mol_o_MO = h_frac * mol_o
    
    # add compounds delivered. This gives total moles of each element in the mantle.
    # from impactor core
    mol_FeO_MO += delivered_fe_s
    mol_ni_MO += delivered_ni
    # from impactor mantle
    mol_si_MO += delivered_si
    mol_v_MO += delivered_v
    mol_mg_MO += delivered_mg
    mol_fe_MO = mol_FeO_MO + delivered_fe
    mol_o_MO += delivered_o
    # Mg and O aren't partitioned into the core, so simply add moles delivered to total moles in the mantle.
    mol_mg += delivered_mg
    mol_o += delivered_o
    # equilibrium and mass balance
    P_eq = Peq(g_acc, h)
    pressure.append(P_eq)
    # print("pressure", P_eq)
    T_eq = Teq(P_eq * 1e9)
    temperature.append(T_eq)
    # print("temperature", T_eq)
    # small error due to exclusion of delta h from volume of impactor mantle (does not change results)
    fe_metal = bisection_search(f, root_bracket(f, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, v_metal), mol_fe_MO, 10e-12, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, v_metal)
    fe_sil = mol_fe_MO - fe_metal
    # both K_eq are decreasing
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.06, 1553, -98)
    ni_sil = mol_ni_MO * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni_MO - ni_sil
    si_sil = (mol_o_MO - ni_sil - fe_sil - mol_mg_MO) / 2
    si_metal = mol_si_MO - si_sil
    v_metal = bisection_search(g, 0, root_bracket(g, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, fe_metal), 10e-12, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, fe_metal)
    v_sil = mol_v_MO - v_metal

    # add moles in metal phase to total moles of each element in the core
    mols_fe_c += fe_metal
    mols_ni_c += ni_metal
    mols_si_c += si_metal
    mols_v_c += v_metal

    # update total moles in mantle of elements involved in equilibrium based on number of moles left in the silicate phase
    mol_fe += fe_sil - mol_fe * h_frac
    mol_ni += ni_sil - mol_ni * h_frac
    mol_si += si_sil - mol_si * h_frac
    mol_v += v_sil - mol_v * h_frac

    # change impactor composition
    fe_s -= 0.02
    si_s = (54 - 0.00606 - fe_s / 2) / (molar_mass_si + molar_mass_o * 2) * molar_mass_si
    mg_s = (46 - fe_s / 2) / (molar_mass_mg + molar_mass_o) * molar_mass_mg
    o_s = 100 - fe_s - mg_s - si_s - v_s

    # add moles in silicate phase to total moles of each element in the mantle
    new_mantle_mass = convert_moles_to_mass(mol_fe, molar_mass_fe) + convert_moles_to_mass(mol_ni, molar_mass_ni) + convert_moles_to_mass(mol_si, molar_mass_si) + convert_moles_to_mass(mol_v, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg)
    new_core_mass = convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_v_c, molar_mass_v)

    # increase planet size
    planet_core_radius = sphere_radius(new_core_mass, rho_core)
    new_mantle_depth = shell_width(new_mantle_mass, rho_mantle, planet_core_radius)
    mantle_depth.append(new_mantle_depth)
    r_planet = new_mantle_depth + planet_core_radius

    # calculate fugacity
    X_FeO.append(mol_fe / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    X_Fe.append(mols_fe_c / (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c))
    fO2.append(calculate_ln_o_iw_fugacity(X_FeO[-1], X_Fe[-1]))

    # track other values
    X_Si.append(mols_si_c / (mols_fe_c + mols_v_c + mols_ni_c + mols_si_c))
    X_NiO.append(mol_ni / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    X_Ni.append(mols_ni_c / (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c))
    X_SiO2.append(mol_si / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    X_Va.append(mols_v_c / (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c))
    X_VO.append(mol_v * 2 / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    X_Mg.append(mol_mg / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    j += 1
save_data(X_Fe, X_Si, X_Ni, X_Va, X_FeO, X_SiO2, X_NiO, X_Mg, X_VO, X_FeO_impactor, gravity, pressure, temperature, planet_size, impactor_size, mantle_depth, fO2, "./data/rapid-solidification/heterogeneous_accretion_decreasing_FeO.txt")
