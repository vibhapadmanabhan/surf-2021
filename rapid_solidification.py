from equilibrium import *
from growth import *
from atmosredox import H2O_H2ratio, GH2O, fO2_fromIW
from atmosphere import Keq_FeO_H2O, GFeO
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

# r_impactor = 10e3
# c_impactor = core_radius(r_impactor)
# delivered_fe = 1000 * fe * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_fe
# delivered_ni = 1000 * ni * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_ni
# delivered_fe_s = fe_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) * 1000 / molar_mass_fe
# delivered_si = 1000 * si_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_si
# delivered_mg = 1000 * mg_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_mg
# delivered_v = 1000 * v_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_v
# delivered_o = 1000 * o_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_o
v_metal = 0
j = 0

for i in range(10):
    print(j)
    planet_size.append(r_planet)
    # if j != 0 and j % 1000 == 0:
    r_impactor = r_planet / 5
    c_impactor = core_radius(r_impactor)
    delivered_fe = 1000 * fe * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_fe
    delivered_ni = 1000 * ni * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_ni
    delivered_fe_s = fe_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) * 1000 / molar_mass_fe
    delivered_si = 1000 * si_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_si
    delivered_mg = 1000 * mg_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_mg
    delivered_v = 1000 * v_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_v
    delivered_o = 1000 * o_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_o
    impactor_size.append(r_impactor)
    X_FeO_impactor.append(fe_s)
    # print("total mantle depth", r_planet - planet_core_radius)
    h = calculate_h(melt_factor * calculate_vol(r_impactor), r_planet)
    # print("melted depth", h)
    # h_frac is the fraction of the mantle that is molten (magma ocean)
    h_frac = calculate_shell_vol(h, r_planet) / calculate_shell_vol(r_planet - planet_core_radius, r_planet)
    planet_mass = convert_moles_to_mass(mol_fe, molar_mass_fe) + convert_moles_to_mass(mol_ni, molar_mass_ni) + convert_moles_to_mass(mol_si, molar_mass_si) + convert_moles_to_mass(mol_v, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg) + convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_v_c, molar_mass_v)
    g_acc = calculate_g(planet_mass)
    gravity.append(g_acc)

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
    T_eq = Teq(P_eq * 1e9)
    temperature.append(T_eq)
    # small error due to exclusion of delta h from volume of impactor mantle (does not change results)
    fe_metal = bisection_search("fe", root_bracket("fe", mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, v_metal), mol_fe_MO, 10e-12, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, v_metal)
    fe_sil = mol_fe_MO - fe_metal
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.06, 1553, -98)
    actual_kd_si = calculate_kd("si", T_eq, P_eq, 2.98, -15934, 0)
    ni_sil = mol_ni_MO * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni_MO - ni_sil
    si_sil = (mol_o_MO - ni_sil - fe_sil - mol_mg_MO) / 2
    si_metal = mol_si_MO - si_sil
    v_metal = bisection_search("v", 0, root_bracket("v", mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, fe_metal), 10e-12, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, fe_metal)
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

    new_mantle_mass = convert_moles_to_mass(mol_fe, molar_mass_fe) + convert_moles_to_mass(mol_ni, molar_mass_ni) + convert_moles_to_mass(mol_si, molar_mass_si) + convert_moles_to_mass(mol_v, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg)
    new_core_mass = convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_v_c, molar_mass_v)

    # increase planet size
    planet_core_radius = sphere_radius(new_core_mass, rho_core)
    new_mantle_depth = shell_width(new_mantle_mass, rho_mantle, planet_core_radius)
    mantle_depth.append(new_mantle_depth)
    r_planet = new_mantle_depth + planet_core_radius

    # calculate fugacity
    X_FeO.append(mol_fe / (mol_fe + mol_ni + mol_si + mol_v / 2 + mol_mg))
    # X_Fe.append(mols_fe_c / (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c))
    X_Fe.append(fe_metal / (ni_metal + si_metal + fe_metal + v_metal))
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

# atmosphere
# convert fO2 to bars
fO2_bar = fO2_fromIW(np.exp(fO2[-1]), Teq(0))
print("fO2 in bars", fO2_bar)

# calculate moles of H2O assuming H2O is 0.1 wt% of mantle weight and all H2 is in the form of water
mol_H2O = 0.01 * h_s * new_mantle_mass * 1000 / (molar_mass_h * 2 + molar_mass_o)
mol_h = h_frac * (2 * mol_H2O)

# total moles of Fe available for re-equilibrium is number of moles left in MO after impact
mol_fe_reeq = fe_sil

# total moles of O
mol_o_atmos = fe_sil + h_frac * mol_H2O

# ignoring CO2 for now
# mol_CO2 = 0.001 * new_mantle_mass * 1000 / (molar_mass_c + molar_mass_o * 2)

# # equilibrium reaction 2H2 + O2 <--> 2H2O occurs
# # using fO2 from equilibrium calculate H2O / H2 ratio
# H2O_H2 = H2O_H2ratio(fO2_bar, Teq(0))

# # equilibrium reaction FeO + H2 <--> Fe + H2O occurs
# # calculate concentration of Fe in MO using equilibrium relation
# xFeO_MO = H2O_H2 / Keq_FeO_H2O(Teq(0))

# fe sil value is 6.103461777373827e+17
mol_fe_mo = bisection_search("mol_fe_mo", root_bracket("mol_fe_mo", fe_sil, ni_sil, si_sil, mol_o_atmos, v_sil, mol_mg * h_frac, P_eq, T_eq, mol_h), fe_sil, 1e-9, fe_sil, ni_sil, si_sil, mol_o_atmos, v_sil, mol_mg * h_frac, 0, Teq(0), mol_h)
mol_fe_metal = mol_fe - mol_fe_mo
mol_H2O = mol_o - mol_fe_mo
mol_H2 = (mol_h - 2 * mol_H2O) / 2
conc_fe_mo = mol_fe_mo / (ni_sil + mol_fe_mo + si_sil + v_sil + h_frac * mol_mg)
conc_H2O = mol_H2O / (mol_H2O + mol_H2)
conc_H2 = mol_H2 / (mol_H2O + mol_H2)

print("calculated kd", conc_H2O / conc_H2 / conc_fe_mo)
conc_fe_metal = 1
final_fO2 = calculate_ln_o_iw_fugacity(conc_fe_mo, conc_fe_metal)
print(final_fO2)

# print(fO2[-10:])
# save_data(X_Fe, X_Si, X_Ni, X_Va, X_FeO, X_SiO2, X_NiO, X_Mg, X_VO, X_FeO_impactor, gravity, pressure, temperature, planet_size, impactor_size, mantle_depth, fO2, "./data/rapid-solidification/.txt")
