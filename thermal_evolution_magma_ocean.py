import matplotlib.pyplot as plt
import math

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
    print(F_conv)
    dT_m = -F_conv * dt / c_e_per_unit_mass / rho_0 / (r_earth - L)
    T_m += dT_m
    L = 1450 * T_m - 2900000
    t += dt

plt.plot(times, temps)
plt.xlabel("Time (yr)")
plt.ylabel("T_m (K)")
plt.title("MO temperature vs time")
plt.show()
plt.plot(times, depths)
plt.xlabel("Time (yr)")
plt.ylabel("Mantle Depth (m)")
plt.title("Mantle depth vs time")
plt.show()
plt.plot(depths, temps)
plt.show()
