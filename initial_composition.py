# starting mantle composition of oxygen-free elements, and oxygen, in wt% (assuming FeO is 8 wt%, MgO is 42 wt%, SiO2 is 50 - 0.00606 wt%, V is 0.00606 wt%). From Lyubetskaya and Korenaga, then scaled to 100. About 8% of mass was missing without scaling. 1 : 1.0868
scale_factor = 1.0866
fe_s = 6.22 * scale_factor # FeO wt% = 8.00
mg_s = 23.41 * scale_factor# 38.82
si_s = 21.09 * scale_factor# 45.19
ni_s = 0
v_s = 0.00606 * scale_factor

# from Hirschmann 2009 assuming most water rich starting composition
h_s = 0.01
c_s = 0.01

# rest is oxygen
o_s = 100 - fe_s - mg_s - si_s - ni_s - v_s - h_s - c_s


# starting core composition
fe = 85
ni = 15
si = 0
v = 0

# molar masses (g / mol)
molar_mass_fe = 55.8
molar_mass_ni = 58.7
molar_mass_si = 28
molar_mass_o = 16
molar_mass_v = 50.9415
molar_mass_mg = 24.3
molar_mass_h = 1
molar_mass_c = 12