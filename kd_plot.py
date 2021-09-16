import csv
import pandas as pd
import matplotlib.pyplot as plt
from equilibrium import *
from growth import *

f = "./data/rapid-solidification/melt_factor_10.txt"
df = pd.read_csv(f, sep='\t', header=0)
kd_si = []
for P in df["P (GPa)"]:
    kd_si.append(calculate_kd("si", Teq(P), P, 2.98, -15934, 0))

plt.plot(df["P (GPa)"], kd_si)
plt.xlabel("Pressure (GPa)")
plt.ylabel("K_d(Si-Fe)")
plt.title("K_d(Si-Fe) dependence on pressure (rapid solidification)")
plt.show()

f = "./data/deep-MO/homogeneous_accretion.txt"
df = pd.read_csv(f, sep='\t', header=0)
kd_si = []
for P in df["P (GPa)"]:
    kd_si.append(calculate_kd("si", Teq(P), P, 2.98, -15934, 0))

plt.plot(df["P (GPa)"], kd_si)
plt.xlabel("Pressure (GPa)")
plt.ylabel("K_d(Si-Fe)")
plt.title("K_d(Si-Fe) dependence on pressure (deep MO)")
plt.show()