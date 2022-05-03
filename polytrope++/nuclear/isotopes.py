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
    str = str.tolower()
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

def AZ_to_str(AZ):
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
        print(value, AZ, AZ == value)
        if AZ == value:
            return key

    A, Z = AZ
    return elements[Z - 1] + str(A)

#def reaction(reaction):
#    string = f"{reaction['reactant'][0]}({reaction['reactant'][1:]},{reaction['product'][1:]}){reaction['product'][0]}"
#    print(string)
    #{"reactant":(), "product"}
