from fe_equilibrium import *
from growth import *

planet_core_radius = core_radius(r_planet)
# in the planet
mol_fe = fe_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) * 1000 / (molar_mass_fe)
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
i = 0
while r_planet <= 1000e3:
    print(i)
    planet_size.append(r_planet)
    r_impactor = r_planet / 10
    c_impactor = core_radius(r_impactor)
    delivered_fe = 1000 * fe * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_fe
    delivered_fe_s = fe_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) * 1000 / (molar_mass_fe)
    delivered_ni = 1000 * ni * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_ni
    delivered_si = 1000 * si_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / (molar_mass_si)
    delivered_mg = 1000 * mg_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / (molar_mass_mg)
    delivered_v = 1000 * v_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / (molar_mass_v)
    delivered_o = 1000 * o_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / (molar_mass_o)
    h = calculate_h(melt_factor * calculate_vol(r_impactor), r_planet)
    print(h)
    mantle_depth.append(h)
    # h_frac is the fraction of the mantle that is molten (magma ocean)
    h_frac = calculate_shell_vol(h, r_planet) / calculate_shell_vol(r_planet - planet_core_radius, r_planet)
    planet_mass = convert_moles_to_mass(mol_fe, molar_mass_fe) + convert_moles_to_mass(mol_ni, molar_mass_ni) + convert_moles_to_mass(mol_si, molar_mass_si) + convert_moles_to_mass(mol_v, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg) + convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_v_c, molar_mass_v)
    g_acc = calculate_g(planet_mass)
    gravity.append(g_acc)
    
    print("gravity", g_acc)
    
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
    print("pressure", P_eq)
    T_eq = Teq(P_eq * 1e9)
    temperature.append(T_eq)
    print("temperature", T_eq)
    # small error due to exclusion of delta h from volume of impactor mantle (does not change results)
    fe_metal = bisection_search(f, root_bracket(f, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, v_metal), mol_fe_MO, 10e-7, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, v_metal)
    fe_sil = mol_fe_MO - fe_metal
    # both K_eq are decreasing
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.06, 1553, -98)
    # actual_kd_si = calculate_kd("si", T_eq, P_eq, 2.98, -15934, None)
    ni_sil = mol_ni_MO * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni_MO - ni_sil
    si_sil = (mol_o_MO - ni_sil - fe_sil - mol_mg_MO) / 2
    si_metal = mol_si_MO - si_sil
    # print("all values", fe_metal, fe_sil, ni_metal, ni_sil, si_metal, si_sil)
    v_metal = bisection_search(g, 0, root_bracket(g, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, fe_metal), 10e-7, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, fe_metal)
    v_sil = mol_v_MO - v_metal

    # add moles in metal phase to total moles of each element in the core
    mols_fe_c += fe_metal
    mols_ni_c += ni_metal
    mols_si_c += si_metal
    mols_v_c += v_metal
    # print("si added molar ratio", si_metal / (fe_metal + ni_metal + si_metal + v_metal))

    # add moles in silicate phase to total moles of each element in the mantle
    added_fe_mantle = fe_sil - h_frac * mol_fe
    added_ni_mantle = ni_sil - h_frac * mol_ni
    added_v_mantle = v_sil - h_frac * mol_v
    added_si_mantle = si_sil - h_frac * mol_si

    # calculate added masses (to core and mantle)
    added_mantle_mass = convert_moles_to_mass(added_fe_mantle, molar_mass_fe) + convert_moles_to_mass(added_ni_mantle, molar_mass_ni) + convert_moles_to_mass(added_si_mantle, molar_mass_si) + convert_moles_to_mass(added_v_mantle, molar_mass_v) + convert_moles_to_mass(delivered_o, molar_mass_o) + convert_moles_to_mass(delivered_mg, molar_mass_mg)
    added_core_mass = convert_moles_to_mass(fe_metal, molar_mass_fe) + convert_moles_to_mass(ni_metal, molar_mass_ni) + convert_moles_to_mass(si_metal, molar_mass_si) + convert_moles_to_mass(v_metal, molar_mass_v)
    added_core_depth = added_depth(planet_core_radius, added_core_mass, rho_core)

    # increase planet size
    planet_core_radius += added_core_depth
    r_planet += added_core_depth
    added_mantle_depth = added_depth(r_planet, added_mantle_mass, rho_mantle)
    r_planet += added_mantle_depth
    print(r_planet)

    # update total moles in mantle of elements involved in equilibrium based on number of moles left in the silicate phase
    mol_fe += fe_sil - mol_fe * h_frac
    mol_ni += ni_sil - mol_ni * h_frac
    mol_si += si_sil - mol_si * h_frac
    mol_v += v_sil - mol_v * h_frac

    # calculate fugacity
    X_FeO.append(mol_fe / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    X_Fe.append(mols_fe_c / (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c))
    print("molar conc of Fe in core", X_Fe)
    fO2.append(calculate_ln_o_iw_fugacity(X_FeO[-1], X_Fe[-1]))

    # track other values
    X_Si.append(mols_si_c / (mols_fe_c + mols_v_c + mols_ni_c + mols_si_c))
    X_NiO.append(mol_ni / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    X_Ni.append(mols_ni_c / (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c))
    X_SiO2.append(mol_si / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    X_Va.append(mols_v_c / (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c))
    X_VO.append(mol_v * 2 / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    i += 1

save_data(X_FeO, X_Fe, X_SiO2, X_Si, X_Va, X_VO, X_Ni, X_NiO, pressure, temperature, gravity, planet_size, fO2, mantle_depth)

plt.plot([r / 1e3 for r in planet_size], fO2)
plt.xlabel("Planet radius (km)")
plt.ylabel("ln(fO2) - IW")
plt.title("Fugacity vs. planet radius")
plt.show()
