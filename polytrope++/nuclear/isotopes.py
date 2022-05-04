"""Utility functions for workig with isotopes and nuclear reactions
"""


elements = ["h", "he",
            "li", "be", "b", "c", "n", "o", "f", "ne",
            "na", "mg", "al", "si", "p", "s", "cl", "ar",
            "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni"]

particles = {"n":(1, 0), "p":(1, 1), "d":(2, 1), "t":(3, 1), "a":(4, 2)}

def str_to_AZ(str):
    """isotope string to (A, Z)

    Take string of the form element abreviation, mass number or particle
    abreviation and convert it to (A, Z).

    Parameters
    ----------
    str : string
        Isotope or particle string, ex "fe56" or "p"

    Returns
    -------
    tuple
        The tuple (A, Z) where A is the mass number and Z is the atomic
        number.
    """
    str = str.lower()
    if str in particles:
        return particles[str]

    A, Z = -1, -1
    for i in range(1, len(str)):
       if str[i:].isnumeric():
           A = int(str[i:])
           if str[:i] in elements:
                Z = elements.index(str[:i]) + 1
                break
    return (A, Z)

def AZ_to_str(AZ, capitalize=True):
    """(A, Z) to isotope string

    (A, Z) to string of the form element abreviation, mass number or
    particle abreviation.

    Parameters
    ----------
    AZ : tuple
        Tuple (A, Z) where A is the mass number and Z is the atomic
        number.

    Returns
    -------
    string
        Isotope or particle string, ex "fe56" or "p"
    """
    for key, value in particles.items():
        if AZ == value:
            return key

    A, Z = AZ
    if capitalize:
        return elements[Z - 1].capitalize() + str(A)
    else:
        return elements[Z - 1] + str(A)

def reaction_name(reaction):
    """Get reaction name

    Parameters
    ----------
    reaction : dict
        dict of reactants and products

    Returns
    -------
    string
        Reaction equation
    """
    reactants = [AZ_to_str(reactant) for reactant in reaction["reactant"]]
    products = [AZ_to_str(product) for product in reaction["product"]]
    return f"{reactants[0]}({' '.join(reactants[1:])},{' '.join(products[1:])}){products[0]}"
