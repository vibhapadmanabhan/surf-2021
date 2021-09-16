from equilibrium import *
from growth import *


"""This file models the deep magma ocean concept. The depth of the MO is always 70% the depth of the mantle."""

planet_core_radius = core_radius(r_planet)


# molar amounts of elements in planetesimal's mantle
mol_fe = 1000 * fe_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_fe)
mol_ni = 1000 * ni_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_ni)
mol_si = 1000 * si_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_si)
mol_mg = 1000 * mg_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_mg)
mol_v = 1000 * v_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_v)
mol_o = 1000 * o_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_o)


# molar amounts of elements in planetesimal's core
mols_si_c = 1000 * si * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_si
mols_v_c = 1000 * v * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_v
mols_ni_c = 1000 * ni * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_ni
mols_fe_c = 1000 * fe * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_fe


# set to 0.
v_metal = 0
j = 0


while r_planet <= 1000e3: 
    print(j)
    planet_size.append(r_planet)
    X_FeO_impactor.append(fe_s)

    # calculate delivered amounts of material
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

    # h_frac is the fraction of the mantle that is molten (magma ocean). 
    h_frac = calculate_shell_vol(0.7 * (r_planet - planet_core_radius), r_planet) / calculate_shell_vol(r_planet - planet_core_radius, r_planet)

    # calculate the mass of the planet
    planet_mass = convert_moles_to_mass(mol_fe, molar_mass_fe) + convert_moles_to_mass(mol_ni, molar_mass_ni) + convert_moles_to_mass(mol_si, molar_mass_si) + convert_moles_to_mass(mol_v, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg) + convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_v_c, molar_mass_v)

    # and its gravitational acceleration
    g_acc = calculate_g(planet_mass)
    gravity.append(g_acc)
    

    # calculate compounds already present in mantle up till melted depth assuming a homogeneous mantle
    # mantle_mass = convert_moles_to_mass(mol_fe, molar_mass_fe) + convert_moles_to_mass(mol_ni, molar_mass_ni) + convert_moles_to_mass(mol_si, molar_mass_si) + convert_moles_to_mass(mol_v, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg)
    h = 0.7 * (r_planet - planet_core_radius)


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


    # calculate current pressure and temperature
    P_eq = Peq(g_acc, h)
    # pressure.append(P_eq)
    T_eq = Teq(P_eq * 1e9)
    # temperature.append(T_eq)


    # Solve using equilibrium and mass balance equations for the amount of iron in metal phase.
    fe_metal = bisection_search("fe", root_bracket("fe", mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, v_metal), mol_fe_MO, 1e-2, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, v_metal)
    fe_sil = mol_fe_MO - fe_metal
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.06, 1553, -98)
    actual_kd_si = calculate_kd("si", T_eq, P_eq, 2.98, -15934, 0)
    ni_sil = mol_ni_MO * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni_MO - ni_sil
    si_sil = (mol_o_MO - ni_sil - fe_sil - mol_mg_MO) / 2
    si_metal = mol_si_MO - si_sil
    v_metal = bisection_search("v", 0, root_bracket('v', mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, fe_metal), 10e-12, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, fe_metal)
    v_sil = mol_v_MO - v_metal

    mantle_solid_depth = r_planet - planet_core_radius - h
    metal_pond_depth = shell_width(fe_metal, 8000, planet_core_radius + mantle_solid_depth)
    r1 = planet_core_radius + mantle_solid_depth
    r2 = r1 + metal_pond_depth
    print("metal pond thickness in km", metal_pond_depth / 1e3)
    pond_timescale = metal_pond_timescale(r1, r2, mantle_solid_depth, planet_core_radius)

    metal_pond_lifetime.append(pond_timescale)
    print("pond timescale", pond_timescale)
    # # add moles in metal phase to total moles of each element in the melt pond (which eventually sinks down into the planetary core)
    mols_fe_c += fe_metal
    mols_ni_c += ni_metal
    mols_si_c += si_metal
    mols_v_c += v_metal

    # update total moles in mantle of elements involved in equilibrium based on number of moles left in the silicate phase
    mol_fe += fe_sil - mol_fe * h_frac
    mol_ni += ni_sil - mol_ni * h_frac
    mol_si += si_sil - mol_si * h_frac
    mol_v += v_sil - mol_v * h_frac


    # calculate new mantle and core mass
    new_mantle_mass = convert_moles_to_mass(fe_sil, molar_mass_fe) + convert_moles_to_mass(ni_sil, molar_mass_ni) + convert_moles_to_mass(si_sil, molar_mass_si) + convert_moles_to_mass(v_sil, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg)
    new_core_mass = convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_v_c, molar_mass_v)

    # increase planet size
    planet_core_radius = sphere_radius(new_core_mass, rho_core)
    new_mantle_depth = shell_width(new_mantle_mass, rho_mantle, planet_core_radius)
    mantle_depth.append(new_mantle_depth)
    r_planet = planet_core_radius + new_mantle_depth


    # calculate ln fO2 - IW
    X_FeO.append(fe_sil / (fe_sil + ni_sil + si_sil + v_sil / 2 + mol_mg * h_frac))
    # X_Fe.append(mols_fe_c / (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c))
    X_Fe.append(fe_metal / (ni_metal + si_metal + fe_metal + v_metal))
    fO2.append(calculate_ln_o_iw_fugacity(X_FeO[-1], X_Fe[-1]))
    # print(fO2[-1])


    # track other values
    # X_Si.append(mols_si_c / (mols_fe_c + mols_v_c + mols_ni_c + mols_si_c))
    X_Si.append(si_metal / (si_metal + fe_metal + ni_metal + v_metal))
    X_NiO.append(ni_sil / (fe_sil + ni_sil + si_sil + mol_v / 2 + mol_mg * h_frac))
    X_Ni.append(ni_metal / (fe_metal + si_metal + ni_metal + v_metal))
    X_SiO2.append(si_sil / (fe_sil + ni_sil + si_sil + v_sil / 2 + mol_mg * h_frac))
    X_Va.append(v_metal/ (fe_metal + si_metal + ni_metal + v_metal))
    X_VO.append(v_sil * 2 / (fe_sil + ni_sil + si_sil + v_sil / 2 + mol_mg * h_frac))
    X_Mg.append(mol_mg * h_frac / (fe_sil + ni_sil + si_sil + v_sil / 2 + mol_mg * h_frac))
    j += 1


feo_wt_percents = [i * (molar_mass_fe + molar_mass_o) * (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg)/ 1000 / new_mantle_mass for i in X_FeO[-10:]]
print(feo_wt_percents)

si_wt_percents = [i * (molar_mass_si) * (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c)/ 1000 / new_core_mass for i in X_Si[-10:]]
print(si_wt_percents)
# save_data(X_Fe, X_Si, X_Ni, X_Va, X_FeO, X_SiO2, X_NiO, X_Mg, X_VO, X_FeO_impactor, gravity, pressure, temperature, planet_size, impactor_size, mantle_depth, fO2, "./data/deep-MO/homogeneous_accretion_earth_size_onefifth_r_planet_impactors_70percent_CMBpressure_fullmantleequilibration.txt")
