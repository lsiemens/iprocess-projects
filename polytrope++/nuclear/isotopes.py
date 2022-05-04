"""Utility functions for workig with isotopes and nuclear reactions
"""


elements = ["h", "he",
            "li", "be", "b", "c", "n", "o", "f", "ne",
            "na", "mg", "al", "si", "p", "s", "cl", "ar",
            "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni"]

particles = {"n":(1, 0), "p":(1, 1), "d":(2, 1), "t":(3, 1), "a":(4, 2)}

def str_to_AZN(str):
    """isotope string to (A, Z, N)

    Take string of the form element abreviation, mass number or particle
    abreviation and convert it to (A, Z, N).

    Parameters
    ----------
    str : string
        Isotope or particle string, ex "fe56" or "2p"

    Returns
    -------
    tuple
        The tuple (A, Z, N) where A is the mass number, Z is the atomic
        number and N is the number of particles.
    """
    n = len(str)
    N = 1
    for i in range(n):
        # shift index by one to check from n-1 to 0 chars at the start
        # of the string
        if str[:n - i - 1].isnumeric():
            N = int(str[:n - i - 1])
            str = str[n - i - 1:]
            break

    str = str.lower()
    if str in particles:
        return (*particles[str], N)

    A, Z = -1, -1
    for i in range(1, len(str)):
       if str[i:].isnumeric():
           A = int(str[i:])
           if str[:i] in elements:
                Z = elements.index(str[:i]) + 1
                break
    if (A < 0) or (Z < 0):
        raise RuntimeError("Could not get mass number or atomic number")
    return (A, Z, N)

def AZN_to_str(AZN, capitalize=True):
    """(A, Z, N) to isotope string

    (A, Z, N) to string of the form element abreviation, mass number or
    particle abreviation.

    Parameters
    ----------
    AZ : tuple
        Tuple (A, Z, N) where A is the mass number, Z is the atomic
        number and N is the number of particles.

    Returns
    -------
    string
        Isotope or particle string, ex "fe56" or "p"
    """
    A, Z, N = AZN
    if N == 1:
        num = ""
    else:
        num = str(N)

    for key, value in particles.items():
        if (A, Z) == value:
            return num + key

    elem = elements[Z - 1]
    if capitalize:
        elem = elem.capitalize()
    mass = str(A)
    return num + elem + mass

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
    # get reactants if more than one particel of main reactant split it
    # into two and move the first reactant to main_reactant, then
    # convert AZN to string
    reactants = reaction["reactant"]
    if reactants[0][-1] != 1:
        main_reactant = (*reactants[0][:-1], 1)
        A, Z, N = reactants[0]
        reactants[0] = (A, Z, N - 1)
    else:
        main_reactant = reactants[0]
        reactants = reactants[1:]
    main_reactant = AZN_to_str(main_reactant)
    reactants = " ".join([AZN_to_str(reactant) for reactant in reactants])

    products = reaction["product"]
    if products[0][-1] != 1:
        main_product = (*products[0][:-1], 1)
        A, Z, N = products[0]
        products[0] = (A, Z, N - 1)
    else:
        main_product = products[0]
        products = products[1:]
    main_product = AZN_to_str(main_product)
    products = " ".join([AZN_to_str(product) for product in products])

    return main_reactant + "(" + reactants + ", " + products + ")" + main_product
