from equilibrium import *
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
for i in range(10):
    print(i)
    planet_size.append(r_planet)
    r_impactor = r_planet / 10
    c_impactor = core_radius(r_impactor)
    delivered_fe = 1000 * fe * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_fe
    delivered_ni = 1000 * ni * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_ni
    delivered_fe_s = fe_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) * 1000 / molar_mass_fe
    delivered_si = 1000 * si_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_si
    delivered_mg = 1000 * mg_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_mg
    delivered_v = 1000 * v_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_v
    delivered_o = 1000 * o_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_o
    h = r_planet - planet_core_radius

    planet_mass = convert_moles_to_mass(mol_fe, molar_mass_fe) + convert_moles_to_mass(mol_ni, molar_mass_ni) + convert_moles_to_mass(mol_si, molar_mass_si) + convert_moles_to_mass(mol_v, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg) + convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_v_c, molar_mass_v)
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
    fe_metal = bisection_search(f, root_bracket(f, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, v_metal), mol_fe, 10e-7, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, v_metal)
    fe_sil = mol_fe - fe_metal
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.06, 1553, -98)
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil - mol_mg) / 2
    si_metal = mol_si - si_sil
    v_metal = bisection_search(g, 0, root_bracket(g, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, fe_metal), 10e-7, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, fe_metal)
    v_sil = mol_v - v_metal

    # add moles in metal phase to total moles of each element in the core
    mols_fe_c += fe_metal
    mols_ni_c += ni_metal
    mols_si_c += si_metal
    mols_v_c += v_metal

    # add moles in silicate phase to total moles of each element in the mantle
    added_fe_mantle = fe_sil - mol_fe
    added_ni_mantle = ni_sil - mol_ni
    added_v_mantle = v_sil - mol_v
    added_si_mantle = si_sil - mol_si

    # calculate added masses (to core and mantle)
    added_mantle_mass = convert_moles_to_mass(added_fe_mantle, molar_mass_fe) + convert_moles_to_mass(added_ni_mantle, molar_mass_ni) + convert_moles_to_mass(added_si_mantle, molar_mass_si) + convert_moles_to_mass(added_v_mantle, molar_mass_v) + convert_moles_to_mass(delivered_o, molar_mass_o) + convert_moles_to_mass(delivered_mg, molar_mass_mg)
    added_core_mass = convert_moles_to_mass(fe_metal, molar_mass_fe) + convert_moles_to_mass(ni_metal, molar_mass_ni) + convert_moles_to_mass(si_metal, molar_mass_si) + convert_moles_to_mass(v_metal, molar_mass_v)
    added_core_depth = added_depth(planet_core_radius, added_core_mass, rho_core)

    # increase planet size
    planet_core_radius += added_core_depth
    r_planet += added_core_depth
    added_mantle_depth = added_depth(r_planet, added_mantle_mass, rho_mantle)
    r_planet += added_mantle_depth

    # update total moles in mantle of elements involved in equilibrium based on number of moles left in the silicate phase
    mol_fe += added_fe_mantle
    mol_ni += added_ni_mantle
    mol_si += added_si_mantle
    mol_v += added_v_mantle

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
    i += 1

save_deep_MO_data(X_FeO, X_Fe, X_SiO2, X_Si, X_Va, X_VO, X_Ni, X_NiO, pressure, temperature, gravity, planet_size, fO2, mantle_depth)

# plot fugacity
plt.plot([r / 1e3 for r in planet_size], fO2)
plt.xlabel("Planet radius (km)")
plt.ylabel("ln(fO2) - IW")
plt.title("Fugacity vs. planet radius")
plt.show()

# plot Si molar concentration in core
plt.plot([r / 1e3 for r in planet_size], X_Si)
plt.xlabel("Planet radius (km)")
plt.ylabel("[Si]_m")
plt.title("Molar concentration of Si in core vs. planet radius")
plt.show()

# plot FeO molar concentration in mantle
plt.plot([r / 1e3 for r in planet_size], X_FeO)
plt.xlabel("Planet radius (km)")
plt.ylabel("[FeO]_sil")
plt.title("Molar concentration of FeO in mantle vs. planet radius")
plt.show()

