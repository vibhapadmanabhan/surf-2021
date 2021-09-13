from equilibrium import *
from growth import *
from atmosredox import H2O_H2ratio, GH2O, fO2_fromIW
from atmosphere import CO2_COratio, Keq_FeO_H2O, GFeO

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
    # assume all H delivered is in the form of H2O
    delivered_h = 2 * 1000 * h_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / (molar_mass_h * 2 + molar_mass_o)
    # assume all C delivered is in the form of CO2
    delivered_c = 1000 * c_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / (molar_mass_c + molar_mass_o * 2)
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
    X_FeO.append(fe_sil / (fe_sil + ni_sil + si_sil + v_sil / 2 + mol_mg * h_frac))
    # X_Fe.append(mols_fe_c / (mols_fe_c + mols_si_c + mols_ni_c + mols_v_c))
    X_Fe.append(fe_metal / (ni_metal + si_metal + fe_metal + v_metal))
    fO2.append(calculate_ln_o_iw_fugacity(X_FeO[-1], X_Fe[-1]))

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

    # atmosphere
    # convert fO2 to bars
    fO2_bar = fO2_fromIW(np.exp(fO2[-1]), Teq(0))

    # assume CO2 and H2O always make up the same wt % of the MO + impactor system because volatiles delivered and lost at a constant rate
    # calculate moles of H in the magma ocean, assuming all H is present in the form of H2O.
    mol_h = 2 * h_frac * 0.01 * h_s * new_mantle_mass * 1000 / (molar_mass_o + molar_mass_h * 2) + delivered_h
    # calculate moles of C in the magma ocean, assuming all C is present in the form of CO2 (assumption subject to change).
    mol_c = h_frac * 0.01 * c_s * new_mantle_mass * 1000 / (molar_mass_c + molar_mass_o * 2) + delivered_c

    # equilibrium reaction 2H2 + O2 <--> 2H2O occurs
    # using fO2 from equilibrium calculate H2O / H2 ratio
    H2O_H2 = H2O_H2ratio(fO2_bar, Teq(0))
    mol_H2 = 1 / (1 + H2O_H2) * mol_h / 2
    mol_H2O = (mol_h - (2 * mol_H2)) / 2
    # total moles of Fe available for re-equilibrium is number of moles left in MO after impact
    mol_fe_reeq = fe_sil
    # total moles of O available for re-equilibrium
    mol_o_atmos = mol_fe_reeq + mol_H2O
    # # equilibrium reaction FeO + H2 <--> Fe + H2O occurs
    val0 = (mol_o_atmos - mol_h / 2) * 1.000000001
    mol_fe_mo = bisection_search("mol_fe_mo", val0, fe_sil, 1e-3, fe_sil, ni_sil, si_sil, mol_o_atmos, v_sil, mol_mg * h_frac, 0, Teq(0), mol_h)
    mol_fe_metal = mol_fe_reeq - mol_fe_mo
    mol_H2O = mol_o_atmos - mol_fe_mo
    mol_H2 = (mol_h - 2 * mol_H2O) / 2
    conc_fe_mo = mol_fe_mo / (ni_sil + mol_fe_mo + si_sil + v_sil + h_frac * mol_mg)
    conc_H2O = mol_H2O / (mol_H2O + mol_H2)
    conc_H2 = mol_H2 / (mol_H2O + mol_H2)

    # add metal to core and remove from mantle
    
    mol_fe_reeq -= mol_fe_metal

    fO2_reeq = calculate_fugacity(conc_H2O, conc_H2, Keq_FeO_H2O(Teq(0)))
    # print("ln fO2 H2O_H2", math.log(fO2_reeq))

    # equilibrium reaction 2CO + O2 <--> 2CO2
    CO2_CO = CO2_COratio(fO2_reeq, Teq(0))
    mol_CO = 1 / (1 + CO2_CO) * mol_c
    mol_CO2 = mol_c - mol_CO
    mol_o_atmos = mol_fe_mo + mol_CO + mol_CO2 * 2

    # FeO + CO <--> Fe + CO2
    mol_fe_mo = bisection_search("mol_c_mo", root_bracket("mol_c_mo", mol_fe_reeq, ni_sil, si_sil, mol_o_atmos, v_sil, mol_mg * h_frac, 0, Teq(0), mol_c), mol_fe_reeq, 1e-3, mol_fe_reeq, ni_sil, si_sil, mol_o_atmos, v_sil, mol_mg * h_frac, 0, Teq(0), mol_c)
    
    mol_fe_metal += mol_fe_reeq - mol_fe_mo
    mol_carbon_compounds = mol_o - mol_fe_mo
    conc_fe_mo = mol_fe_mo / (mol_ni + mol_fe_mo + mol_si + mol_v + mol_mg)
    CO2_CO = Keq_FeO_CO2(Teq(0)) * conc_fe_mo
    mol_CO = mol_carbon_compounds * 1 / (1 + CO2_CO)
    mol_CO2 = mol_carbon_compounds - mol_CO
    conc_CO = mol_CO / (mol_CO + mol_CO2)
    conc_CO2 = mol_CO2 / (mol_CO + mol_CO2)

    fO2_reeq = calculate_fugacity(conc_CO2, conc_CO, Keq_FeO_CO2(Teq(0)))
    # print("ln fO2 CO2_CO", fO2_reeq)

    # allow the Fe metal to go to core and update size of planet
    new_core_mass = convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_v_c, molar_mass_v)
    planet_core_radius = sphere_radius(new_core_mass, rho_core)

    
    mols_fe_c += mol_fe_metal
    # update Fe in mantle
    mol_fe -= mol_fe_metal
