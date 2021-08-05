from equilibrium import *
"""Generates a plot of fO2-IW w.r.t impactor radius. 
Mantle composition: 10.7 wt% Si, 0 wt% Fe, 0% Mg, 0% V."""
r_impactor = [i * 10 * 1e3 for i in range(1, 101)]
X_FeO = []
X_Fe = []
X_Si = []
ln_fO2 = []

for i in range(len(r_impactor)):
    mol_si = calculate_total_element_moles(si, molar_mass_si, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
    mol_fe = calculate_total_element_moles(fe, molar_mass_fe, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
    mol_ni = calculate_total_element_moles(ni, molar_mass_ni, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
    mol_o = 2 * mol_si
    fe_metal = bisection_search(root_bracket(mol_fe, mol_ni, mol_si, mol_o, 0, 0, P_eq, T_eq, 0), mol_fe, 10e-7, mol_fe, mol_ni, mol_si, mol_o, 0, 0, P_eq, T_eq, 0)
    fe_sil = mol_fe - fe_metal
    actual_kd_ni = calculate_kd('ni', T_eq, P_eq, 1.06, 1553, -98)
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil) / 2
    si_metal = mol_si - si_sil
    X_Si.append(si_metal)
    X_FeO.append(fe_sil / (fe_sil + ni_sil + si_sil))
    X_Fe.append(fe_metal / (fe_metal + ni_metal + si_metal))

r_impactor = [i / 1e3 for i in r_impactor]
plt.plot(r_impactor, X_FeO)
plt.xlabel("Radius of impactor (km)")
plt.ylabel("X_FeO")
plt.title("X_FeO vs impactor radius (for a 1000 km radius planet)")
plt.show()

for i in range(len(r_impactor)):
    ln_fO2.append(calculate_ln_o_iw_fugacity(X_FeO[i], X_Fe[i]))

plt.plot(r_impactor, ln_fO2)
plt.xlabel("Radius of impactor (km)")
plt.ylabel("ln(fO2)")
plt.title("ln(fO2) vs impactor radius (for a 1000 km radius planet)")
plt.show()