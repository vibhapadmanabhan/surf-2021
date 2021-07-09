import numpy as np
# values calculated in fe_equilibrium.py
mol_fe = 7.656916027693821e+19
mol_ni = 1.284465119787499e+19
mol_si = 3.1590935698224997e+22
mol_o = 6.318187139644999e+22

actual_kd_si = 0.21798434487972793
kd_ni = 10.459236904631686

# search for roots of equation
a = 1e17 # negative value for difference between calculated and actual K_d(Si-Fe)
b = mol_fe # positive

eps = 1e-7

def f(fe_metal):
    """
    Returns the difference between the actual value of K_d(Si-Fe) and calculated value for K_d(Si-Fe). (Root of function will be found at the 
    correct value of K_d)
    """
    fe_sil = mol_fe - fe_metal
    ni_sil = mol_ni * fe_sil / (fe_sil + kd_ni * fe_metal)
    ni_metal = mol_ni - ni_sil
    si_sil = (mol_o - ni_sil - fe_sil) / 2
    si_metal = mol_si - si_sil
    return si_metal * fe_sil**2 / si_sil / fe_metal**2 - actual_kd_si

while True:
    FA = f(a)
    fe_metal = (a + b) / 2
    FP = f(fe_metal)
    if np.abs(FP) <= eps: # close enough to true value
        break
    if FA * FP > 0:
        a = fe_metal
    else:
        b = fe_metal


fe_sil = mol_fe - fe_metal
ni_sil = mol_ni * fe_sil / (fe_sil + kd_ni * fe_metal)
ni_metal = mol_ni - ni_sil
si_sil = (mol_o - ni_sil - fe_sil) / 2
si_metal = mol_si - si_sil
print(fe_metal / (fe_metal + ni_metal + si_metal))