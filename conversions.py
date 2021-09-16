
def convert_mass_to_moles(compound_mass, compound_molar_mass):
    """Converts mass to moles given the mass of a compound in g and its molar mass in g / mol."""
    return compound_mass / compound_molar_mass


def convert_moles_to_mass(mol_element, element_molar_mass):
    """Returns the mass in kg of an element given its molar mass and number of moles."""
    return mol_element * element_molar_mass / 1000