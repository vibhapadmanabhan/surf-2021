from growth import *

"""This file models the rapid solidification theory. The magma ocean solidifies completely before the next impactor arrives."""


planet_core_radius = core_radius(r_planet)
melt_factor = 10  # melt volume produced upon impact = melt_factor * vol_impactor
h = r_planet - planet_core_radius  # mantle depth


# molar amounts of elements in planetesimal's mantle
mol_fe = 1000 * fe_s * 0.01 * ocean_mass(r_planet, h) / molar_mass_fe
mol_ni = 1000 * ni_s * 0.01 * ocean_mass(r_planet, h) / molar_mass_ni
mol_si = 1000 * si_s * 0.01 * ocean_mass(r_planet, h) / molar_mass_si
mol_mg = 1000 * mg_s * 0.01 * ocean_mass(r_planet, h) / molar_mass_mg
mol_v = 1000 * v_s * 0.01 * ocean_mass(r_planet, h) / molar_mass_v
mol_o = 1000 * o_s * 0.01 * ocean_mass(r_planet, h) / molar_mass_o

# molar amounts of elements in planetesimal's core
mols_si_c = 1000 * si * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_si
mols_v_c = 1000 * v * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_v
mols_ni_c = 1000 * ni * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_ni
mols_fe_c = 1000 * fe * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_fe


# molar amounts of volatiles in planetesimal's atmosphere
mol_H2O_tot = 0
mol_H2_tot = 0
mol_CO2_tot = 0
mol_CO_tot = 0
mol_CH4_tot = 0


# set to 0
v_metal = 0
j = 0


while r_planet <= 1000e3:  # growing the planetesimal to a planetary embryo
    print(j)
    planet_size.append(r_planet)

    # ========================   MAGMA OCEAN   =================================
    # calculate delivered amounts of material
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
    X_FeO_impactor.append(fe_s) # in case of heterogeneous accretion, track FeO

    # calculate the depth h of the magma ocean formed
    h = calculate_h(melt_factor * calculate_vol(r_impactor), r_planet)

    # h_frac is the fraction of the mantle that is molten (magma ocean)
    h_frac = calculate_shell_vol(h, r_planet) / calculate_shell_vol(r_planet - planet_core_radius, r_planet)


    # calculate the mass of the planet
    planet_mass = convert_moles_to_mass(mol_fe, molar_mass_fe) + convert_moles_to_mass(mol_ni, molar_mass_ni) + convert_moles_to_mass(mol_si, molar_mass_si) + convert_moles_to_mass(mol_v, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg) + convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_v_c, molar_mass_v)

    # and its gravitational acceleration
    g_acc = calculate_g(planet_mass)
    gravity.append(g_acc)

    # calculate compounds already present in mantle up till melted depth assuming a homogeneous mantle
    mol_si_MO = h_frac * mol_si
    mol_v_MO = h_frac * mol_v
    mol_FeO_MO = h_frac * mol_fe
    mol_ni_MO = h_frac * mol_ni
    mol_mg_MO = h_frac * mol_mg
    mol_o_MO = h_frac * mol_o
    
    # add the molar amounts of compounds delivered. This gives total moles of each element in the magma ocean.
    # from impactor core
    mol_FeO_MO += delivered_fe_s
    mol_ni_MO += delivered_ni
    # from impactor mantle
    mol_si_MO += delivered_si
    mol_v_MO += delivered_v
    mol_mg_MO += delivered_mg
    mol_fe_MO = mol_FeO_MO + delivered_fe
    mol_o_MO += delivered_o

    # Mg and O aren't partitioned into the metal, so simply add moles delivered to total moles in the mantle.
    mol_mg += delivered_mg
    mol_o += delivered_o

    # calculate current pressure and temperature
    P_eq = Peq(g_acc, h)
    pressure.append(P_eq)
    T_eq = Teq(P_eq * 1e9)  # pressure needs to be in Pa
    temperature.append(T_eq)


    # Solve using equilibrium and mass balance equations for the amount of iron in metal phase.
    fe_metal = bisection_search("fe", root_bracket("fe", mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, v_metal), mol_fe_MO, 1e-6, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, v_metal)
    fe_sil = mol_fe_MO - fe_metal
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.06, 1553, -98)
    actual_kd_si = calculate_kd("si", T_eq, P_eq, 2.98, -15934, 0)
    ni_sil = mol_ni_MO * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni_MO - ni_sil
    si_sil = (mol_o_MO - ni_sil - fe_sil - mol_mg_MO) / 2
    si_metal = mol_si_MO - si_sil
    print("si metal", si_metal)
    print("fe metal conc", fe_metal / (fe_metal + si_metal + v_metal + ni_metal))
    v_metal = bisection_search("v", 0, root_bracket("v", mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, fe_metal), 1e-6, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, fe_metal)
    v_sil = mol_v_MO - v_metal
    print("completed initial equilibrium")

    # add moles in metal phase to total moles of each element in the melt pond (which eventually sinks down into the planetary core)
    mols_fe_c += fe_metal
    mols_ni_c += ni_metal
    mols_si_c += si_metal
    mols_v_c += v_metal

    # update total moles in mantle of elements involved in equilibrium based on number of moles left in the silicate phase
    mol_fe += fe_sil - mol_fe * h_frac
    mol_ni += ni_sil - mol_ni * h_frac
    mol_si += si_sil - mol_si * h_frac
    mol_v += v_sil - mol_v * h_frac


    # recalculate core and mantle masses
    new_mantle_mass = convert_moles_to_mass(mol_fe, molar_mass_fe) + convert_moles_to_mass(mol_ni, molar_mass_ni) + convert_moles_to_mass(mol_si, molar_mass_si) + convert_moles_to_mass(mol_v, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg)
    new_core_mass = convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_v_c, molar_mass_v)

    # increase planet size
    planet_core_radius = sphere_radius(new_core_mass, rho_core)
    new_mantle_depth = shell_width(new_mantle_mass, rho_mantle, planet_core_radius)
    mantle_depth.append(new_mantle_depth)
    r_planet = new_mantle_depth + planet_core_radius


    # calculate ln(fO2 - IW)
    X_FeO.append(fe_sil / (fe_sil + ni_sil + si_sil + v_sil / 2 + mol_mg * h_frac))
    X_Fe.append(fe_metal / (ni_metal + si_metal + fe_metal + v_metal))
    fO2.append(calculate_ln_o_iw_fugacity(X_FeO[-1], X_Fe[-1]))
    print(fO2[-1])


    # track other values
    X_Si.append(si_metal / (si_metal + fe_metal + ni_metal + v_metal))
    X_NiO.append(ni_sil / (fe_sil + ni_sil + si_sil + mol_v / 2 + mol_mg * h_frac))
    X_Ni.append(ni_metal / (fe_metal + si_metal + ni_metal + v_metal))
    X_SiO2.append(si_sil / (fe_sil + ni_sil + si_sil + v_sil / 2 + mol_mg * h_frac))
    X_Va.append(v_metal/ (fe_metal + si_metal + ni_metal + v_metal))
    X_VO.append(v_sil * 2 / (fe_sil + ni_sil + si_sil + v_sil / 2 + mol_mg * h_frac))
    X_Mg.append(mol_mg * h_frac / (fe_sil + ni_sil + si_sil + v_sil / 2 + mol_mg * h_frac))
    j += 1


    # ========================   ATMOSPHERE   =================================
    # convert fO2 to bars, assuming temperature of the system is the temperature at 0 GPa pressure.
    fO2_bar = fO2_fromIW(np.exp(fO2[-1]), Teq(0))

    # calculate moles of H in the magma ocean, assuming all H is present in the form of H2O.
    mol_h = 2 * h_frac * 0.01 * h_s * new_mantle_mass * 1000 / (molar_mass_o + molar_mass_h * 2) + delivered_h

    # calculate moles of C in the magma ocean, assuming all C is present in the form of CO2.
    mol_c = h_frac * 0.01 * c_s * new_mantle_mass * 1000 / (molar_mass_c + molar_mass_o * 2) + delivered_c

    # add moles of element already in atmosphere
    mol_h += mol_CH4_tot * 4 + mol_H2O_tot * 2 + mol_H2_tot * 2
    mol_c += mol_CH4_tot + mol_CO2_tot + mol_CO_tot


    # using fO2 from equilibrium calculate initial molar amounts of H2O, H2, CO2, CO, CH4 to determine how much oxygen is available for re-equilibrium
    mol_o_atmos = mol_c * 2 + mol_h / 2 + mol_CO2_tot * 2 + mol_CO_tot + mol_H2O_tot

    # set an upper bound for moles of water such that no other amounts become negative
    upper_bound = H2O_H2ratio(fO2_bar, Teq(0)) * mol_h / 2 / (1 + H2O_H2ratio(fO2_bar, Teq(0)))


    # solve for initial amounts of each volatile
    mol_H2O = bisection_search_atmosphere(1e-16, upper_bound, 1e-3, mol_o_atmos, mol_h, mol_c, fO2_bar, Teq(0))
    print("completed oxidation")
    mol_H2 = mol_H2O / H2O_H2ratio(fO2_bar, Teq(0))
    mol_CH4 = (mol_h - 2 * mol_H2 - 2 * mol_H2O) / 4
    mol_CO2 = mol_o_atmos - mol_c - mol_H2O + mol_CH4 # from C and O mass bal.
    mol_CO = mol_CO2 / CO2_COratio(fO2_bar, Teq(0))
    mol_volatiles = mol_H2O + mol_H2 + mol_CO + mol_CO2 + mol_CH4
    mol_fe_reeq = fe_sil
    mol_CH4_tot = mol_CH4  # mol_CH4 is now the total amount of CH4 in atmos.


    # calculate total amount of oxygen in atmosphere-MO system
    mol_o_atmos = mol_fe_reeq + mol_H2O + mol_CO + mol_CO2 * 2
    mol_h = mol_H2 * 2 + mol_H2O * 2   # ignore CH4


    # re-equilibrium between FeO and volatiles
    mol_fe_mo = bisection_search("mol_fe_mo", 1e-6, fe_sil, 1e-6, fe_sil, ni_sil, si_sil, mol_o_atmos, v_sil, mol_mg * h_frac, 0, Teq(0), mol_h)
    kd_CO2 = Keq_FeO_CO2(Teq(0))
    print("keq co2", kd_CO2)
    # kd_CH4 = Keq_FeO_CH4(T_eq)
    kd_H2O = Keq_FeO_H2O(Teq(0))
    print("keq H2o", kd_H2O)
    mol_fe_metal = mol_fe - mol_fe_mo
    conc_fe_mo = mol_fe_mo / (ni_sil + mol_fe_mo + si_sil + mol_mg * h_frac + v_sil)
    H2O_to_H2 = kd_H2O * conc_fe_mo
    mol_H2 = 1 / (H2O_to_H2 + 1) * mol_h / 2 # ignore CH4 here because its concentration is low
    mol_H2O = (mol_h - 2 * mol_H2) / 2
    CO2_to_CO = kd_CO2 * conc_fe_mo
    mol_CO = 1 / (2 * CO2_to_CO + 1) * (mol_o_atmos - mol_fe_mo - mol_H2O)
    mol_CO2 = mol_CO * CO2_to_CO
    # mol_CH4 = mol_c - mol_CO - mol_CO2
    mol_volatiles = mol_H2O + mol_H2 + mol_CO + mol_CO2
    # conc_CH4 = mol_CH4 / mol_volatiles
    conc_CO = mol_CO / mol_volatiles
    print("conc CO2 after FeO re-eq", conc_CO)
    conc_CO2 = mol_CO2 / mol_volatiles
    conc_H2O = mol_H2O / mol_volatiles
    conc_H2 = mol_H2 / mol_volatiles

    mol_CO_tot = mol_CO
    mol_CO2_tot = mol_CO2
    mol_H2_tot = mol_H2
    mol_H2O_tot = mol_H2O

    # allow the Fe metal to go to planetary core and update size of planet
    mols_fe_c += mol_fe_metal
    # remove metal phase Fe from mantle
    mol_fe -= mol_fe_metal


    # once again calculate core mass, mantle mass, and planet radius
    new_core_mass = convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_v_c, molar_mass_v)
    ew_mantle_mass = convert_moles_to_mass(mol_fe, molar_mass_fe) + convert_moles_to_mass(mol_ni, molar_mass_ni) + convert_moles_to_mass(mol_si, molar_mass_si) + convert_moles_to_mass(mol_v, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg)

    # increase planet size
    planet_core_radius = sphere_radius(new_core_mass, rho_core)
    new_mantle_depth = shell_width(new_mantle_mass, rho_mantle, planet_core_radius)
    mantle_depth.append(new_mantle_depth)
    r_planet = new_mantle_depth + planet_core_radius
