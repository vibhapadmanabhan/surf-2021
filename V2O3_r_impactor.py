from fe_equilibrium import *

r_impactor = [i * 10 * 1e3 for i in range(1, 100)]
# r_impactor = 100 * 1e3
P_cmb = 0.7 * 25 # GPa
P_eq = [i * 0.1 * P_cmb for i in range(1, 10)]
X_FeO = []
X_Fe = []
X_Si = []
X_Va = []
v_metal = 0
T_eq = 4000
# mantle_depth = [i * 1e3 for i in range(300, 701)]
for i in range(len(r_impactor)):
    mol_si = calculate_total_element_moles(si, molar_mass_si, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
    mol_fe_s = calculate_total_element_moles(fe_s, molar_mass_fe + molar_mass_o, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
    mol_fe = calculate_total_element_moles(fe, molar_mass_fe, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i])) + mol_fe_s
    mol_ni = calculate_total_element_moles(ni, molar_mass_ni, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
    # consider vanadium in the mantle of the planet
    mol_v = calculate_total_element_moles(v, molar_mass_v, planet_mantle_depth, r_planet, r_impactor[i], calculate_impactor_core_radius(r_impactor[i]))
    mol_o = 2 * mol_si + mol_fe_s # ignore contributions of V oxides in mass balance for O

    fe_metal = bisection_search(f, root_bracket(f, mol_fe, mol_ni, mol_si, mol_o, mol_v, P_eq, T_eq, v_metal), mol_fe, 10e-7, mol_fe, mol_ni, mol_si, mol_o, mol_v, P_eq, T_eq, v_metal)
    actual_kd_ni = calculate_kd("ni", T_eq, P_eq, 1.06, 1553, -98)
    fe_sil = mol_fe - fe_metal
    ni_sil = mol_ni * fe_sil / (fe_sil + actual_kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil) / 2
    si_metal = mol_si - si_sil

    v_metal = bisection_search(g, 0, root_bracket(g, mol_fe, mol_ni, mol_si, mol_o, mol_v, P_eq, fe_metal), 10e-7, mol_fe, mol_ni, mol_si, mol_o, mol_v, P_eq, fe_metal)
    v_sil = mol_v - v_metal

    X_Si.append(si_metal / (fe_metal + ni_metal + si_metal + v_metal))
    X_FeO.append(fe_sil / (fe_sil + ni_sil + si_sil + v_sil))
    X_Fe.append(fe_metal / (fe_metal + ni_metal + si_metal + v_metal))
    X_Va.append(v_metal / (v_metal + ni_metal + fe_metal + si_metal))

r_impactor = [i / 1e3 for i in r_impactor]

fO2 = []
for i in range(len(r_impactor)):
    fO2.append(calculate_ln_o_iw_fugacity(X_FeO[i], X_Fe[i]))

plt.plot(r_impactor, X_Si)
plt.xlabel("Impactor radius (km)")
plt.ylabel("X_V2O3")
plt.title("X_V2O3 vs impactor radius (for a 1000 km radius planet)")
plt.show()
