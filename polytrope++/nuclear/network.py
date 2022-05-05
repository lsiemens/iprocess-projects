"""Manage nuclear reaction networks
"""

import isotopes
from matplotlib import pyplot

def load_network(fname):
    """Load network file

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

    print(text)
    return [isotopes.str_to_reaction(line) for line in text]

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
            print(isotopes.AZN_to_str((*orig, 1)))
            name = isotopes.AZN_to_str((*orig, 1))
            pyplot.gca().annotate(f"${name}$", xy=final, xycoords="data", size=16)
    pyplot.xlabel("$Z$ (element)")
    pyplot.ylabel("$A$ (isotope)")
    pyplot.title("Nuclear reaction network")
    pyplot.gca().set_aspect(1)
    pyplot.grid()
    pyplot.show()

network = load_network("./networks/cnocycle") + load_network("./networks/ppchain")

draw_network(network, True, subs=[[(1, 1), (1, 10)],[(4, 2), (10,1)]])
