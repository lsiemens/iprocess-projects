"""Plots and graphs
"""

from matplotlib import pyplot

from . import isotopes
from . import network

def network_graph(reactions, show_all=False, subs=[], show=True, ifig=None):
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
        If True, show the network. The default is False
    ifig : integer, optional
        Figure number. The default is None
    """
    if ifig is not None:
        pyplot.close(ifig)
        pyplot.figure(ifig)

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
    if show:
        if ifig is not None:
            pyplot.show(ifig)
        else:
            pyplot.show()

def isotope_evolution(t, Y, particles, show=True, ifig=None):
    """Plot evolution of mass fractions

    Parameters
    ----------
    t : array
        The time in [s]
    Y : array
        Molar abundance evolution array in [mol/g]
    particles : list
        List of AZN tuples of the particles in the reaction network
    show : bool, optional
        If True, show the network. The default is False
    ifig : integer, optional
        Figure number. The default is None
    """

    t = t/(60*60*24*365.25)
    X = network.molar_abu_to_mass_frac(particles, Y)

    if ifig is not None:
        pyplot.close(ifig)
        pyplot.figure(ifig)

    for X_i, AZ_i in zip(X, particles):
        label = isotopes.AZN_to_str((*AZ_i, 1), particle_alt=False)
        pyplot.loglog(t, X_i, label=label)
    pyplot.xlabel("t [yr]")
    pyplot.ylabel("mass fraction")
    pyplot.title("Nuclear reaction network evolution")
    pyplot.ylim(1e-10, 1.5)
    pyplot.legend()

    if show:
        if ifig is not None:
            pyplot.show(ifig)
        else:
            pyplot.show()
