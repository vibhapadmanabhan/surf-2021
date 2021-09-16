import math
from planetesimal_initial_parameters import k
from densities import *
import deprecation

def calculate_g(M_p):
    """Calculate gravitational acceleration using g ~ M_p^(0.503) and M_Mars = 6.39e23 kg, g_Mars = 3.7 m / s"""
    return M_p**0.503 * k


def calculate_h(melt_vol, r_planet):
    """Returns the depth of mantle that is melted by an impactor."""
    return r_planet - (r_planet**3 - 3 / 4 / math.pi * melt_vol)**(1 / 3)


def Peq(g, h):
    """Returns P_eq in GPa."""
    return rho_mantle * g * h / 1e9


def added_depth(initial_radius, added_mass, rho):
    """Return increase in radius / depth of core / mantle respectively given added mass."""
    return (3 / 4 / math.pi * added_mass / rho + initial_radius**3)**(1 / 3) - initial_radius


def core_radius(body_radius):
    """Returns initial core radius upon the assumption that initial metal mass is 25% of the total mass."""
    return (body_radius**3 / 9)**(1 / 3)


@deprecated
def core_radius2(body_radius):
    """Returns initial core radius upon the assumption that initial metal mass is 20% of the total mass."""
    return (body_radius**3 / 9.4)**(1 / 3)
    
    
def Teq(Peq):
    """Input Peq in Pa."""
    return 0.4 * 1661.2 * (Peq / 1.336 / 10**9 + 1)**(1 / 7.437) + 0.6 * 1982.1 * (Peq / 6.594 / 1e9 + 1)**(1 / 5.374)