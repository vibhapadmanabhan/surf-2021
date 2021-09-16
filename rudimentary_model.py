from equilibrium import *


"""Generates a plot of fO2-IW w.r.t impactor radius. """
r_impactor = [i * 10 * 1e3 for i in range(1, 101)]  # impactor size in m

X_FeO = []
X_Fe = []
X_Si = []
ln_fO2 = []
X_Va = []
v_metal = 0

for i in range(len(r_impactor)):
    # moles of species from mantles
    mol_si = calculate_total_element_moles(si_s, molar_mass_si, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
    mol_mg = calculate_total_element_moles(mg_s, molar_mass_mg, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
    mol_fe_s = calculate_total_element_moles(fe_s, molar_mass_fe, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
    mol_v = calculate_total_element_moles(v_s, molar_mass_v, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))

    # from impactor core
    mol_fe = calculate_total_element_moles(fe, molar_mass_fe, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
    
    mol_fe += mol_fe_s  # total moles of iron
    mol_ni = calculate_total_element_moles(ni, molar_mass_ni, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
   
    mol_o = calculate_total_element_moles(o_s, molar_mass_o, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))

    fe_metal = bisection_search("fe", root_bracket("fe", mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, v_metal), mol_fe, 1e-12, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, v_metal)

    fe_sil = mol_fe - fe_metal
    actual_kd_ni = calculate_kd('ni', T_eq, P_eq, 1.06, 1553, -98)
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil - mol_mg) / 2
    si_metal = mol_si - si_sil
    v_metal = bisection_search("v", 0, root_bracket('v', mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, fe_metal), 10e-12, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, fe_metal)
    v_sil = mol_v - v_metal

    X_Si.append(si_metal / (si_metal + ni_metal + fe_metal + v_metal))
    X_FeO.append(fe_sil / (fe_sil + ni_sil + si_sil + v_sil + mol_mg))
    X_Fe.append(fe_metal / (fe_metal + ni_metal + si_metal + v_metal))
    X_Va.append(v_sil / (fe_sil + ni_sil + si_sil + v_sil))
    ln_fO2.append(calculate_ln_o_iw_fugacity(X_FeO[-1], X_Fe[-1]))


# print(X_Fe)
# r_impactor = [i / 1e3 for i in r_impactor]
# print(r_impactor[:10])
# plt.plot(r_impactor, X_FeO)
# plt.xlabel("Radius of impactor (km)")
# plt.ylabel("X_FeO")
# plt.title("X_FeO vs impactor radius (for a 1000 km radius planet)")
# plt.show()

r_impactor = [i / 1e3 for i in r_impactor]  # convert to km

# make a plot. Code adapted from https://www.py4u.net/discuss/15394
plt.style.use('tex')
fig = plt.figure()
ax1 = fig.add_subplot(2,1,1)
ax1.plot(r_impactor, ln_fO2, color='k', linewidth=0.7)
ax1.set(ylabel="$\ln f\mathrm{O}_2 - \mathrm{IW}$")
ax1.set_title("(a)", loc="right")
ax1.tick_params(direction='in', length=3,   which='major', top=True, right=True)
ax1.tick_params(direction='in', length=1.5, which='minor', top=True, right=True)
ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)
ax2.plot(r_impactor, X_Va, color='k', linewidth=0.7)
ax2.set(ylabel="$X_V^{sil}$")
ax2.set_title("(b)", loc="right")
ax2.tick_params(direction='in', length=3,   which='major', top=True, right=True)
ax2.tick_params(direction='in', length=1.5, which='minor', top=True, right=True)
ax2.plot(r_impactor, X_Va, color='k', linewidth=0.7)
plt.xlabel("Radius of impactor (km)")
plt.setp(ax1.get_xticklabels(), visible=False)
plt.tight_layout()
plt.show()
