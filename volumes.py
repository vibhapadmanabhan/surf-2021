import math

def sphere_radius(mass, density):
    """Calculates the radius of a sphere, given mass and density"""
    return (3 / 4/ math.pi * mass / density)**(1 / 3)

def shell_width(mass, density, inner_radius):
    """Calculates the width of a spherical shell given the inner radius of the shell, the mass and the density"""
    return (3 / 4 / math.pi * mass / density + inner_radius**3)**(1 / 3) - inner_radius

def shell_width_from_outer(mass, density, outer_radius):
    """Calculates the width of a spherical shell given the outer radius of the shell, the mass and the density."""
    return outer_radius - (- 3 / 4 / math.pi * mass / density + outer_radius**3)**(1 / 3)

def calculate_shell_vol(mantle_depth, r_obj):
    """Returns the volume of a spherical shell given the width of the shell and its outer radius."""
    return 4 / 3 * math.pi * (r_obj**3 - (r_obj - mantle_depth)**3)


def calculate_vol(r_obj):
    """Returns the volume of a sphere."""
    return 4 / 3 * math.pi * r_obj**3