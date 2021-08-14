from equilibrium import *
"""Generates a plot of fO2-IW w.r.t amount of FeO in the mantle. No vanadium and magnesium in mantle."""
r_impactor = 100 * 1e3
fe_s = [i for i in range(50)]
X_FeO = []
X_Fe = []
X_Si = []
X_V2O3 = []
mol_mg = 0
mol_v = 0

for i in range(len(fe_s)):
    mol_si = calculate_total_element_moles(si, molar_mass_si, planet_mantle_depth, r_planet, r_impactor, calculate_impactor_core_radius(r_impactor))
    mol_fe_s = calculate_total_element_moles(fe_s[i], molar_mass_fe + molar_mass_o, planet_mantle_depth, r_planet, r_impactor, calculate_impactor_core_radius(r_impactor))
    mol_fe = calculate_total_element_moles(fe, molar_mass_fe, planet_mantle_depth, r_planet, r_impactor, calculate_impactor_core_radius(r_impactor)) + mol_fe_s
    mol_ni = calculate_total_element_moles(ni, molar_mass_ni, planet_mantle_depth, r_planet, r_impactor, calculate_impactor_core_radius(r_impactor))
    mol_o = 2 * mol_si + mol_fe_s
    fe_metal = bisection_search("fe", root_bracket("fe", mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, 0), mol_fe, 10e-7, mol_fe, mol_ni, mol_si, mol_o, mol_v, mol_mg, P_eq, T_eq, 0)
    fe_sil = mol_fe - fe_metal
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.06, 1553, -98)
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil) / 2
    si_metal = mol_si - si_sil
    X_Si.append(si_metal)
    X_FeO.append(fe_sil / (fe_sil + ni_sil + si_sil))
    X_Fe.append(fe_metal / (fe_metal + ni_metal + si_metal))

fO2 = []
for i in range(len(fe_s)):
    fO2.append(calculate_ln_o_iw_fugacity(X_FeO[i], X_Fe[i]))

plt.plot(fe_s, fO2)
plt.xlabel("Initial wt percent of FeO in the mantle")
plt.ylabel("ln(fO2)")
plt.title("ln(fO2) vs FeO amount in mantle")
plt.show()
