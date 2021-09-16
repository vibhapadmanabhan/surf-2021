import math
from densities import *


def body_mass(r_body, core_radius):
    return 4 / 3 * math.pi * (core_radius**3 * rho_core + (r_body - core_radius)**3 * rho_mantle)


def ocean_mass(r_planet, depth):
    return 4 / 3 * math.pi * (r_planet**3 - (r_planet - depth)**3) * rho_mantle