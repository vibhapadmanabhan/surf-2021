from fe_equilibrium import *
from growth import *

# starting mantle composition in wt%
fe_s = 0
mg_s = 46
si_s = 54 - 0.0606
ni_s = 0
v_s = 0.0606

# starting core composition in wt%
fe = 85
ni = 15
si = 0
v = 0
mg = 0

# in the planet
mol_fe = fe_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) * 1000 / (molar_mass_fe + molar_mass_o)
mols_fe_c = 1000 * fe * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_fe
mol_ni = 1000 * ni_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_ni + molar_mass_o)
mols_ni_c = 1000 * ni * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_ni
mol_si = si_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_si + 2 * molar_mass_o)
mols_si_c = 1000 * si * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_si
mol_mg = 1000 * mg_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_mg + molar_mass_o)
mols_mg_c = 1000 * mg * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_mg
mol_v = (1000 * v_s * 0.01 * ocean_mass(r_planet, r_planet - planet_core_radius) / (molar_mass_v * 2 + molar_mass_o * 3)) * 2
mols_v_c = 1000 * v * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_v

# moles of stuff in the impactor
delivered_fe = 1000 * fe * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_fe
delivered_fe_s = fe_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) * 1000 / (molar_mass_fe + molar_mass_o)
delivered_ni = 1000 * ni * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_ni
delivered_si = 1000 * si_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / (molar_mass_si + 2 * molar_mass_o)
delivered_mg = 1000 * mg_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / (molar_mass_mg + molar_mass_o)
delivered_v = 1000 * v_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / (molar_mass_v * 2 + 3 * molar_mass_o)

for i in range(10):
    print(i)
    planet_size.append(r_planet)
    h = calculate_h(melt_factor * calculate_vol(r_impactor), r_planet)
    h_frac = calculate_shell_vol(h, r_planet) / calculate_shell_vol(r_planet - planet_core_radius, r_planet)
    # h_frac is the fraction of the mantle that is molten (magma ocean)
    g_acc = calculate_g(body_mass(r_planet, planet_core_radius))
    # calculate compounds already present in mantle up till melted depth assuming uniform distribution of elements over mantle volume
    mol_si_MO = h_frac * mol_si
    mol_v_MO = h_frac * mol_v
    mol_FeO_MO = h_frac * mol_fe
    mol_ni_MO = h_frac * mol_ni
    mol_mg_MO = h_frac * mol_mg
    # add compounds delivered. This gives total moles of each element in the mantle.
    # from impactor core
    mol_FeO_MO += delivered_fe_s
    mol_ni_MO += delivered_ni
    # from impactor mantle
    mol_si_MO += delivered_si
    mol_v_MO += delivered_v
    mol_mg_MO += delivered_mg
    mol_fe_MO = mol_FeO_MO + delivered_fe
    mol_o_MO = 2 * mol_si_MO + mol_FeO_MO + mol_ni_MO + mol_mg_MO
    mol_mg += delivered_mg
    # increase h (account for volume of impactor's mantle)
    added_mantle_depth = added_depth(r_planet, body_mass(r_impactor, c_impactor), rho_mantle)
    h += added_mantle_depth
    print("impactor mass", body_mass(r_impactor, c_impactor))
    # equilibrium and mass balance
    P_eq = Peq(g_acc, h) # * 10
    fe_metal = bisection_search(f, root_bracket(f, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, v_metal), mol_fe_MO, 10e-7, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, v_metal)
    fe_sil = mol_fe_MO - fe_metal
    print("sil, metal", fe_sil, fe_metal)
    actual_kd_ni = calculate_kd("ni", T_cmb, P_eq, 1.06, 1553, -98)
    ni_sil = mol_ni_MO * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni_MO - ni_sil
    si_sil = (mol_o_MO - ni_sil - fe_sil - mol_mg_MO) / 2
    si_metal = mol_si_MO - si_sil
    v_metal = bisection_search(g, 0, root_bracket(g, fe_sil, ni_sil, si_sil, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, fe_metal), 10e-7, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, fe_metal)
    v_sil = mol_v_MO - v_metal

    # calculate fugacity
    mols_fe_c += fe_metal
    mols_ni_c += ni_metal
    mols_si_c += si_metal
    mols_v_c += v_metal
    
    # increase planet size
    added_core_mass = convert_moles_to_mass(fe_metal, molar_mass_fe) + convert_moles_to_mass(ni_metal, molar_mass_ni) + convert_moles_to_mass(si_metal, molar_mass_si) + convert_moles_to_mass(v_metal, molar_mass_v)
    added_core_depth = added_depth(planet_core_radius, added_core_mass, rho_core)

    planet_core_radius += added_core_depth
    r_planet += added_mantle_depth + added_core_depth

    mol_fe += fe_sil - fe_metal - mol_fe * h_frac
    print("fe in mantle after equilibrium", mol_fe)
    mol_ni += ni_sil - ni_metal - mol_ni * h_frac
    mol_si += si_sil - si_metal - mol_si * h_frac
    # print(si_metal / mol_si)
    mol_v += v_sil - v_metal - mol_v * h_frac
    # print(mol_fe, mol_ni, mol_si, mol_v, mol_mg)
    X_FeO = mol_fe / (mol_fe + mol_ni + mol_si + mol_v * 2 + mol_mg)
    #print(X_FeO)
    X_Fe = mols_fe_c / (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c)
    #print(X_Fe)
    X_Si.append(mols_si_c / (mols_fe_c + mols_v_c + mols_ni_c + mols_si_c))
    fO2.append(calculate_ln_o_iw_fugacity(X_FeO, X_Fe))
    print("total mass", body_mass(r_planet, planet_core_radius))


plt.plot([r / 1e3 for r in planet_size], fO2)
plt.xlabel("Planet radius (km)")
plt.ylabel("ln(fO2) - IW")
plt.title("Fugacity vs. planet radius")
plt.show()
