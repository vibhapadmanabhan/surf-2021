import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib import rc

# plt.rcParams['figure.figsize']= 4*1.414, 4*1
# plt.rcParams['font.family']='serif'
# plt.rcParams['font.size']     = 11
# plt.rcParams['axes.linewidth']= 0.5
# rc('font',**{'family':'serif','serif':['Computer Modern']})
# rc('text',usetex=True)
# rc('text.latex', preamble=r'\usepackage{sfmath}')

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

# constants
r_earth = 6400000  # m

# convective heat loss
k_c = 2  # W m^-1 K^-1
c_e_per_unit_mass = 1000  # J kg^-1 K^-1
rho_0 = 3300  # kg m^-3
T_s = 300  # K

# rayleigh number
alpha = 2 * 10**(-5)  # K^-1
g = 9.8  # ms^-1
kappa = 10**(-6)  # thermal diffusivity
nu = 10  # Pa s  (viscosity)

# model
T_m = 4000  # K
L = 2900000 # m
dt = 365 * 24 * 3600
t = 0

# lists
temps = []
depths = []
times = []

while L >= 0:
    temps.append(T_m)
    depths.append(L)
    times.append(t / (365 * 24 * 3600))
    Ra = alpha * rho_0 * g * (T_m - T_s) * L**3 / kappa / nu
    F_conv = 0.089 * k_c * (T_m - T_s) / L * Ra**(1 / 3)
    dT_m = -F_conv * dt / c_e_per_unit_mass / rho_0 / (r_earth - L)
    T_m += dT_m
    L = 1450 * T_m - 2900000
    t += dt


# plot surface temperature evolution

# ax.set_xlim([np.min(l_fO2), np.max(l_fO2)])
# ax.set_ylim([0, 1])


# plot evolution of magma ocean depth
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(times, depths, color='k', linewidth=0.7)
ax1.set(xlabel="Time (yr)", ylabel="Magma ocean depth (m)")
ax1.set_title("(a)", loc="right")
ax1.tick_params(direction='in', length=3,   which='major', top=True, right=True)
ax1.tick_params(direction='in', length=1.5, which='minor', top=True, right=True)
ax2.plot(times, temps, color='k', linewidth=0.7)
ax2.set(xlabel="Time (yr)", ylabel="Temperature (K)")
ax2.set_title("(b)", loc="right")
ax2.tick_params(direction='in', length=3,   which='major', top=True, right=True)
ax2.tick_params(direction='in', length=1.5, which='minor', top=True, right=True)
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
