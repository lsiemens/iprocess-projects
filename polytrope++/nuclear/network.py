"""Manage nuclear reaction networks
"""

import numpy
from matplotlib import pyplot
import requests

import isotopes
import reaclib

def read_network_file(fname):
    """Read network file

    Parameters
    ----------
    fname : string
        Path to network file

    Returns
    -------
    list
        List of nuclear reaction dicts
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

    Returns
    -------
    list
        List of reaclib rate tuples of the form,
        (reaction dict, Q-value, reaction_rate function)
    """
    reactions_path = reaclib.read_reaction_list(dir)

    download = None
    reactions_data = []
    for reaction in reactions:
        reaction = isotopes.reaction_to_str(reaction)
        # try to get reaction rate file path from reaction list
        try:
            path = reactions_path[reaction]
        except KeyError:
            print(f"ERROR: no file for the reaction \"{reaction}\" "
                   "found in reaction list. Try downloading the "
                   "missing reaction from JINA Reaclib and/or "
                   "regenerating the reaction list.\n")
            raise

        # try to read reaction rate data
        try:
            reactions_data.append(reaclib.read_file(path))
        except FileNotFoundError:
            if download is None:
                download = input("Download missing Reaclib files? [y/n]:").lower()
            if not download.startswith("y"):
                raise

            print(f"Downloading rate data for {reaction}")
            # try to download missing reaction rate data from JINA Reaclib
            try:
                reaclib.download_reaction(reaction, reactions_path[reaction])
            except requests.HTTPError:
                # try downloading again with out specifying the dataset
                # for the reaction
                reaclib.download_reaction(reaction, reactions_path[reaction],
                                          use_set=False)
            reactions_data.append(reaclib.read_file(path))
    return reactions_data

def build_network(reactions, T9, rho, dir="./reactions"):
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

    Returns
    -------
    list
        List of AZN tuples of the particles in the reaction network
    function
        Function dYdt(t, Y), is the time derivative of the molar
        abundance vector Y. Note that each molar abundance Y[i]
        corrisponds to the particle with the AZN given in particles[i]
    """
    reactions_data = load_network(reactions, dir)
    particles = []
    for reaction in reactions:
        particles += reaction["reactants"] + reaction["products"]
    particles = [particle[:-1] for particle in particles]
    particles = list(set(particles))

    # sort particles by element with isotopes in order of their mass
    sort_secondary = lambda AZN:AZN[0]
    particles = sorted(particles, key=sort_secondary)
    sort_primary = lambda AZN:AZN[1]
    particles = sorted(particles, key=sort_primary)

    mask_As = []
    N_As = []
    mask_Bs = []
    N_Bs = []
    rates = []

    for i, (reaction, q_value, rate) in enumerate(reactions_data):
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

        mask_As.append(mask_A)
        N_As.append(N_A)

        mask_Bs.append(mask_B)
        N_Bs.append(N_B)

        rates.append(rate(T9))

    def dYdt(t, Y):
        """Reaction network differental equation

        dYdt(t, Y), is the time derivative of the molar
        abundance vector Y. Note that each molar abundance Y[i]
        corrisponds to the particle with the AZN given in particles[i]

        Parameters
        ----------
        time : float
            The time
        Y : array
            Molar abundance vector. Note that Y[i] = X[i]/A_i where X[i]
            is the mass fraction and A_i is the molar mass of the isotope

        Returns
        -------
        array
            The differential dYdt
        """
        Y = numpy.asarray(Y)
        dYdt = numpy.zeros(Y.shape)
        for mask_A, N_A, mask_B, N_B, rate in zip(mask_As, N_As, mask_Bs, N_Bs, rates):
            Y_A = Y[mask_A]
            rate_factor = reaclib.rate_factor(Y_A, N_A, rate, rho)

            dYdt[mask_A] -= N_A*rate_factor
            dYdt[mask_B] += N_B*rate_factor
        return dYdt

    return particles, dYdt

def draw_network(reactions, show_all=False, subs=[]):
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