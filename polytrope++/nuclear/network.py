"""Manage nuclear reaction networks
"""

import numpy
from matplotlib import pyplot

import isotopes
import reaclib

def read_network_file(fname):
    """Read network file

    Parameters
    ----------
    fname : string
        Path to network file
    """

    with open(fname, "r") as fin:
        text = fin.read().split("\n")
    text = [line.strip().lower() for line in text]
    text = [line for line in text if line != ""]
    text = [line for line in text if not line.startswith("#")]

    return [isotopes.str_to_reaction(line) for line in text]

def load_network(reactions, dir="./reactions/"):
    """Load network

    Parameters
    ----------
    reactions : list
        List of nuclear reaction dicts.
    dir : string, optional
        Directory to check for reactions. The default is "./reactions/"
    """

    reactions_path = reaclib.read_reaction_list(dir)

    reactions_data = []
    for reaction in reactions:
        reaction = isotopes.reaction_to_str(reaction)
        try:
            path = reactions_path[reaction]
        except KeyError:
            print(f"ERROR: no file for the reaction \"{reaction}\" "
                  f"found in \"{dir}\". Try downloading the missing "
                   "reaction from JINA Reaclib and/or regenerate the "
                   "reaction list.\n")
            raise

        reactions_data.append(reaclib.read_file(path))
    return reactions_data

def build_network(reactions, T9, rho, dir="./reactions"): # TODO
    """Build network at constant T9 and rho

    Parameters
    ----------
    reactions : list
        List of nuclear reaction dicts.
    T9 : float
        The tempurature in Gigakelvin.
    rho : float
        The density in g/cm^3
    dir : string
        Directory to check for reactions. The default is "./reactions/"
    """
    reactions_data = load_network(reactions, dir)
    particles = []
    for reaction in reactions:
        particles += reaction["reactants"] + reaction["products"]
    particles = [particle[:-1] for particle in particles]
    particles = list(set(particles))

    sort_secondary = lambda AZN:AZN[0]
    particles = sorted(particles, key=sort_secondary)
    sort_primary = lambda AZN:AZN[1]
    particles = sorted(particles, key=sort_primary)

    def dYdt(t, Y):
        Y = numpy.asarray(Y)
        dYdt = numpy.zeros(Y.shape)
        for reaction, q_value, llambda in reactions_data:
            reactants = reaction["reactants"]
            products = reaction["products"]

            mask_A = []
            N_A = []
            for AZN in reactants:
                mask_A.append(particles.index(AZN[:-1]))
                N_A.append(AZN[-1])
            N_A = numpy.array(N_A)

            N_B = []
            mask_B = []
            for AZN in products:
                mask_B.append(particles.index(AZN[:-1]))
                N_B.append(AZN[-1])
            N_B = numpy.array(N_B)

            Y_A = Y[mask_A]
            rate_factor = reaclib.rate_factor(Y_A, N_A, llambda, T9, rho)

            dYdt[mask_A] -= N_A*rate_factor
            dYdt[mask_B] += N_B*rate_factor
        return dYdt

    return particles, dYdt

def draw_network(reactions, show_all=False, subs=[], show=True):
    """Draw nuclear reaction network

    Parameters
    ----------
    reactions : list
        List of nuclear reaction dicts.
    show_all : bool, optional
        If true, show secondary isotopes in the reactions. The default
        is False.
    subs : list, optional
        list of pairs of AZNs to substitue, only applies to secondary
        isotopes. The default is [].
    show : bool, optional
        If true the graph will show automaticaly. The defalut is False.
    """
    for reaction in reactions:
        reactants = reaction["reactants"]
        products = reaction["products"]
        main_reactant = reactants[0]
        main_product = products[0]
        A0, Z0 = main_reactant[:-1]
        A1, Z1 = main_product[:-1]
        pyplot.arrow(Z0, A0, Z1 - Z0, A1 - A0, width=0.05, length_includes_head=True)

        if show_all:
            for reactant in reactants:
                A, Z = reactant[:-1]
                for orig, final in subs:
                    if (A, Z) == orig:
                        A, Z = final
                        break
                pyplot.plot([Z, Z0], [A, A0], "r-", alpha=0.5)
            for product in products:
                A, Z = product[:-1]
                for orig, final in subs:
                    if (A, Z) == orig:
                        A, Z = final
                        break
                pyplot.plot([Z1, Z], [A1, A], "g-", alpha=0.5)
    if show_all:
        for orig, final in subs:
            name = isotopes.AZN_to_str((*orig, 1))
            pyplot.gca().annotate(f"${name}$", xy=final, xycoords="data", size=16)
    pyplot.xlabel("$Z$ (element)")
    pyplot.ylabel("$A$ (isotope)")
    pyplot.title("Nuclear reaction network")
    pyplot.gca().set_aspect(1)
    pyplot.grid()
    pyplot.show()

#reactions = read_network_file("./networks/pp_branch_I")
#particles, dydt = build_network(reactions, 1, 1)
#Y_0 = [0.97, 0.1, 0.1, 0.1]
#print(dydt(0.0, Y_0))
