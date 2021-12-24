from growth import *
from initial_composition2 import *
import matplotlib.pyplot as plt

"""This file models the rapid solidification theory. The magma ocean solidifies completely before the next impactor arrives."""
final_fO2_list = []

planet_core_radius = core_radius(r_planet)
melt_factor = 10  # melt volume produced upon impact = melt_factor * vol_impactor
h = r_planet - planet_core_radius  # mantle depth

# ----- [ initial condition ] -----
# molar amounts of elements in planetesimal's mantle
M_MO = ocean_mass(r_planet, h)
mol_fe = fe_s * 0.01 * M_MO / molar_mass_fe
mol_ni = ni_s * 0.01 * M_MO / molar_mass_ni
mol_si = si_s * 0.01 * M_MO / molar_mass_si
mol_mg = mg_s * 0.01 * M_MO / molar_mass_mg
mol_v =  v_s * 0.01 * M_MO  / molar_mass_v
mol_o =  o_s * 0.01 * M_MO  / molar_mass_o

# molar amounts of elements in planetesimal's core
mols_si_c = si * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_si
mols_v_c =   v * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_v
mols_ni_c = ni * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_ni
mols_fe_c = fe * 0.01 * calculate_vol(planet_core_radius) * rho_core / molar_mass_fe


# molar amounts of volatiles in planetesimal's atmosphere
mol_H2O_tot = 0
mol_H2_tot = 0
mol_CO2_tot = 0
mol_CO_tot = 0
mol_CH4_tot = 0


# set to 0
v_metal = 0
j = 0

# save variable
r_H2O  = []
r_nore = []


while r_planet <= 1000e3:  # growing the planetesimal to a planetary embryo
    print("")
    print(j)
    planet_size.append(r_planet)

    # ========================   MAGMA OCEAN   =================================
    # calculate delivered amounts of material
    r_impactor = r_planet / 5
    c_impactor = core_radius(r_impactor)
    delivered_fe =  fe * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_fe
    delivered_ni =  ni * 0.01 * calculate_vol(c_impactor) * rho_core / molar_mass_ni
    delivered_fe_s = fe_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_fe
    delivered_si =   si_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_si
    delivered_mg =   mg_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_mg
    delivered_v =    v_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_v
    delivered_o =    o_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / molar_mass_o
    # assume all H delivered is in the form of H2O
    delivered_h = 2 * h_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / (molar_mass_h * 2 + molar_mass_o)
    # assume all C delivered is in the form of CO2
    delivered_c =     c_s * 0.01 * ocean_mass(r_impactor, r_impactor - c_impactor) / (molar_mass_c + molar_mass_o * 2)


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
    #print("kd error", f(fe_metal, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, v_metal))
    fe_sil = mol_fe_MO - fe_metal
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.06, 1553, -98)
    actual_kd_si = calculate_kd("si", T_eq, P_eq, 2.98, -15934, 0)
    ni_sil = mol_ni_MO * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni_MO - ni_sil
    si_sil = (mol_o_MO - ni_sil - fe_sil - mol_mg_MO) / 2
    si_metal = mol_si_MO - si_sil
    #print("si metal", si_metal)
    print("fe metal conc", fe_metal / (fe_metal + si_metal + v_metal + ni_metal))
    #print("fe sil, ni sil, fe metal, ni metal , si sil, si metal, mol o MO", fe_sil, ni_sil, fe_metal,ni_metal, si_sil, si_metal, mol_o_MO)
    v_metal = bisection_search("v", 0, root_bracket("v", mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, fe_metal), 1e-6, mol_fe_MO, mol_ni_MO, mol_si_MO, mol_o_MO, mol_v_MO, mol_mg_MO, P_eq, T_eq, fe_metal)
    #print("after vanadium")
    v_sil = mol_v_MO - v_metal
    #print("completed initial equilibrium")

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
    mol_h = 2 * h_frac * 0.01 * h_s * new_mantle_mass  / (molar_mass_o + molar_mass_h * 2) + delivered_h
    
    # calculate moles of C in the magma ocean, assuming all C is present in the form of CO2.
    mol_c = h_frac * 0.01 * c_s * new_mantle_mass  / (molar_mass_c + molar_mass_o * 2) + delivered_c
    
    # assume the redox H & C would be reset by the presence of metal pond
    rH2O_pond   = H2O_H2ratio(fO2_bar, Teq(0))
    rCO2_pond   = CO2_COratio(fO2_bar, Teq(0))
    r_nore.append(rH2O_pond)
    
    mol_H2  = 0.5 * mol_h / (rH2O_pond + 1) # ignore CH4 here because its concentration is low
    mol_H2O = 0.5 * (mol_h - 2*mol_H2)
    mol_CO  = mol_c / (rCO2_pond + 1)
    mol_CO2 = mol_c - mol_CO
    mol_CH4 = 0.0 #(mol_h - 2 * mol_H2 - 2 * mol_H2O) / 4
    
    mol_volatiles = mol_H2O + mol_H2 + mol_CO + mol_CO2 # + mol_CH4
    print("fe, ni, si ", fe_sil, "\t", ni_sil, "\t", si_sil)
    print("initial mol_H2, H2O ", mol_H2, mol_H2O, " CO2: ",mol_CO2, "\t r = ", rH2O_pond, " Fe", fe_sil)
    
    # calculate total amount of oxygen in atmosphere-MO system
    mol_fe_reeq = fe_sil
    mol_o_atmos = mol_fe_reeq + mol_H2O + mol_CO + mol_CO2 * 2
    #print("atmospheric oxygen", mol_o_atmos)
    
    # re-equilibrium between FeO and volatiles
    #print("Teq = ", T_eq)
    xmetal = calc_metal(mol_fe_reeq, ni_sil, si_sil, mol_mg*h_frac, 4000, 0) #T_eq, P_eq)
    mol_fe_mo = bisection_search("mol_fe_mo", mol_fe_reeq*.8, mol_fe_reeq*.99999, 1e-4, fe_sil, ni_sil, si_sil, mol_o_atmos, v_sil, mol_mg * h_frac, mol_c, Teq(0), mol_h)
    kd_CO2 = Keq_FeO_CO2(Teq(0))
    kd_H2O = Keq_FeO_H2O(Teq(0))
    #print("keq H2o", kd_H2O)
    #print("keq co2", kd_CO2)
    mol_fe_metal = mol_fe_reeq - mol_fe_mo

    #
    xmetal = .1
    conc_fe_mo = mol_fe_mo / (ni_sil + mol_fe_mo + si_sil + mol_mg * h_frac + v_sil)
    H2O_to_H2  = kd_H2O * conc_fe_mo / xmetal
    mol_H2     = .5 * mol_h / (H2O_to_H2 + 1)   # ignore CH4 here because its concentration is low
    mol_H2O    = .5 * (mol_h - 2 * mol_H2)
    CO2_to_CO  = kd_CO2 * conc_fe_mo / xmetal
    mol_CO     = mol_c / (CO2_to_CO + 1)
    mol_CO2    = mol_c - mol_CO
    mol_volatiles = mol_H2O + mol_H2 + mol_CO + mol_CO2
    print("percentage of iron removed", mol_fe_metal / mol_fe_reeq, "\t", mol_fe_metal, " mol of Fe created, ", mol_H2O, mol_H2)
    print("r after = ", H2O_to_H2)
    
    #print("after mol CO2, mol CO, mol_CH4, mol_H2, mol_H2O, fe in metal, fe_sil", mol_CO2, mol_CO, mol_CH4, mol_H2, mol_H2O, mol_fe_metal, mol_fe_mo)
    
    print("atmospheric composition", mol_CO, mol_CO2, mol_H2, mol_H2O)
    total_CO.append(mol_CO)
    total_CO2.append(mol_CO2)
    total_H2.append(mol_H2)
    total_H2O.append(mol_H2O)
    total_CH4.append(mol_CH4)

    # Fe transport from MO to core by re-equilibrium
    mols_fe_c += mol_fe_metal        # allow the Fe metal to go to planetary core and update size of planet
    mol_fe    -= mol_fe_metal        # remove metal phase Fe from mantle
    
    # once again calculate core mass, mantle mass, and planet radius
    new_core_mass  = convert_moles_to_mass(mols_fe_c, molar_mass_fe) + convert_moles_to_mass(mols_ni_c, molar_mass_ni) + convert_moles_to_mass(mols_si_c, molar_mass_si) + convert_moles_to_mass(mols_v_c, molar_mass_v)
    new_mantle_mass = convert_moles_to_mass(mol_fe, molar_mass_fe) + convert_moles_to_mass(mol_ni, molar_mass_ni) + convert_moles_to_mass(mol_si, molar_mass_si) + convert_moles_to_mass(mol_v, molar_mass_v) + convert_moles_to_mass(mol_o, molar_mass_o) + convert_moles_to_mass(mol_mg, molar_mass_mg)

    # increase planet size
    planet_core_radius = sphere_radius(new_core_mass, rho_core)
    new_mantle_depth = shell_width(new_mantle_mass, rho_mantle, planet_core_radius)
    mantle_depth.append(new_mantle_depth)
    r_planet = new_mantle_depth + planet_core_radius

    mol_volatiles_tot = mol_H2 + mol_H2O + mol_CO2 + mol_CO + mol_CH4
    fraction_H2O = mol_H2O / mol_volatiles_tot
    fraction_H2  = mol_H2  / mol_volatiles_tot
    final_fO2    = calculate_fugacity(fraction_H2O, fraction_H2, Keq_FeO_H2O(Teq(0)))
    final_fO2_list.append(final_fO2)
    print("final fO2", final_fO2)

    r_H2O.append(mol_H2O/mol_H2) #(mol_H2O+mol_H2))
    
# print(total_H2O)

# ----- [plot results] ------
# fig, ax = plt.subplots()
# ax.plot(np.array(planet_size) / 1e3, r_H2O,  color="k") #total_CO)
# ax.plot(np.array(planet_size) / 1e3, r_nore, color="k", linestyle="--") #total_CO)
# ax.set(xlabel="Planetesimal radius [km]", ylabel="H$_2$O/H$_2$ ratio in the atmosphere")
# plt.tight_layout()
# plt.savefig("./total_.pdf")

tex_fonts = {
    # from https://jwalton.info/Embed-Publication-Matplotlib-Latex/
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 10,
    "font.size": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8
}

plt.rcParams.update(tex_fonts)

# plot evolution of magma ocean depth
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot([i / 1e3 for i in planet_size], melt_pond_timescale, color='k', linewidth=0.7)
ax1.set(xlabel="Time (yr)", ylabel="Magma ocean depth (m)")
ax1.set_title("(a)", loc="right")
ax1.tick_params(direction='in', length=3,   which='major', top=True, right=True)
ax1.tick_params(direction='in', length=1.5, which='minor', top=True, right=True)
# ax2.plot(times, temps, color='k', linewidth=0.7)
# ax2.set(xlabel="Time (yr)", ylabel="Temperature (K)")
# ax2.set_title("(b)", loc="right")
# ax2.tick_params(direction='in', length=3,   which='major', top=True, right=True)
# ax2.tick_params(direction='in', length=1.5, which='minor', top=True, right=True)
plt.tight_layout()
# plt.savefig("./plots/magma_ocean_thermal_evolution.pdf", transparent=True)
plt.show()

# plt.plot(times, temps)
# plt.xlabel("Time (yr)")
# plt.ylabel("T_m (K)")
# plt.title("MO temperature vs time")
# plt.show()
# plt.plot(times, depths)
# plt.xlabel("Time (yr)")
# plt.ylabel("Mantle Depth (m)")
# plt.title("Mantle depth vs time")
# plt.show()
# plt.plot(depths, temps)
# plt.show()

