"""Manage nuclear reaction networks
"""

import numpy
import requests
import scipy.integrate
from matplotlib import pyplot

from iniabu import inimf # use mass fractions

from . import isotopes
from . import reaclib
from . import plot

N_Avogadro = 6.02214076e23 # avagadros number in [#/mol]

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

def load_network(reactions, dir="./nuclear/reactions/"):
    """Load network

    Parameters
    ----------
    reactions : list
        List of nuclear reaction dicts.
    dir : string, optional
        Directory to check for reactions. The default is
        "./nuclear/reactions/"

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

def build_network(reactions, dir="./nuclear/reactions/"):
    """Build network

    Parameters
    ----------
    reactions : list
        List of nuclear reaction dicts.
    dir : string
        Directory to check for reactions. The default is
        "./nuclear/reactions/"

    Returns
    -------
    list
        List of AZN tuples of the particles in the reaction network
    function
        function that gets dYdt and epsilon at defined density and
        tempurature
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
    q_values = []

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

        q_values.append(q_value)

    def get_dYdt_epsilon(rho, T):
        """Get dYdt and epsilon given tempurature and density

        Parameters
        ----------
        rho : float
            The density in [g/cm^3]
        T : float
            The tempurature in [K]

        Returns
        -------
        function
            Function dYdt(t, Y), is the time derivative of the molar
            abundance vector Y. Note that each molar abundance Y[i]
            corrisponds to the particle with the AZN given in particles[i]
        function
            Function epsilon(t, Y), is the nuclear reaction power density
        """
        rates = []
        for i, (reaction, q_value, rate) in enumerate(reactions_data):
            rates.append(rate(T))

        def dYdt(t, Y):
            """Reaction network differental equation

            dYdt(t, Y), is the time derivative of the molar
            abundance vector Y. Note that each molar abundance Y[i]
            corrisponds to the particle with the AZN given in particles[i]

            Parameters
            ----------
            t : float
                The time in [s]
            Y : array
                Molar abundance vector in [mol/g]. Note that Y[i] = X[i]/A_i
                where X[i] is the mass fraction and A_i is the molar
                mass of the isotope

            Returns
            -------
            array
                The differential dYdt in [mol/(s g)]
            """
            Y = numpy.asarray(Y)
            dYdt_value = numpy.zeros(Y.shape)
            for mask_A, N_A, mask_B, N_B, rate in zip(mask_As, N_As, mask_Bs, N_Bs, rates):
                Y_A = Y[mask_A]
                rate_factor = reaclib.rate_factor(Y_A, N_A, rate, rho)

                dYdt_value[mask_A] -= N_A*rate_factor
                dYdt_value[mask_B] += N_B*rate_factor
            return dYdt_value

        def epsilon(t, Y):
            """Nuclear energy production rate per unit mass

            epsilon(t, Y), is the rate of thermonuclear energy production
            per unit mass in [MeV/(s g)]

            Parameters
            ----------
            time : float
                The time in [s]
            Y : array
                Molar abundance vector in [mol/g]. Note that Y[i] = X[i]/A_i
                where X[i] is the mass fraction and A_i is the molar
                mass of the isotope

            Returns
            -------
            float
                Energy production rate per unit mass in [MeV/(s g)]
            """
            Y = numpy.asarray(Y)
            epsilon_value = 0
            for mask_A, N_A, mask_B, N_B, rate, q_value in zip(mask_As, N_As, mask_Bs, N_Bs, rates, q_values):
                Y_A = Y[mask_A]
                rate_factor = reaclib.rate_factor(Y_A, N_A, rate, rho)
                epsilon_value += rate_factor*N_Avogadro*q_value
            return epsilon_value
        return dYdt, epsilon
    return particles, get_dYdt_epsilon

def solve_network(rho, T, Y_0, get_dYdt_epsilon, t_range=(0., 1.), **kwargs):
    """Solve nuclear reaction network

    Parameters
    ----------
    rho : float
        The density in [g/cm^3]
    T : float
        The tempurature in [K]
    Y_0 : array
        The inital molar abundance vector in [mol/g]
    get_dYdt_epsilon : function
        Function that gets dYdt and epsilon at defined density and
        tempurature for a given nuclear network
    t_range : tuple, optional
        Time range in [s]. The default is (0., 1.)
    kwargs : dict
        Key word arguments for scipy solve_ivp

    Returns
    -------
    t : array
        time coordinate array in [s]
    Y : array
        Molar abundance evolution array in [mol/g]
    status : integer
        Status of solve_ivp
    """
    dYdt, _ = get_dYdt_epsilon(rho, T)
    solve = scipy.integrate.solve_ivp(dYdt, t_range, Y_0, method="Radau", **kwargs)
    if not solve.success:
        print(solve.message)
    return solve["t"], solve["y"], solve.status

def solve_radial_network(rho, T, Y_0, get_dYdt_epsilon, t_range=(0., 1.), plot_args=None, **kwargs):
    """Solve nuclear reaction network along a radial profile

    Parameters
    ----------
    rho : array
        The density in [g/cm^3]
    T : array
        The tempurature in [K]
    Y_0 : array
        The inital molar abundance vector in [mol/g]
    get_dYdt_epsilon : function
        Function that gets dYdt and epsilon at defined density and
        tempurature for a given nuclear network
    t_range : tuple, optional
        Time range in [s]. The default is (0., 1.)
    kwargs : dict
        Key word arguments for scipy solve_ivp

    Returns
    -------
    Y : array
        Molar abundance evolution array in [mol/g]
    status : integer
        Status of solve_ivp
    """
    if plot_args is not None:
        particles, show, ifig = plot_args

    Y = numpy.zeros(shape=Y_0.shape)
    t = numpy.zeros(shape=rho.shape)
    status = 0
    for i in range(len(rho)):
        t_i, Y_i, status_i = solve_network(rho[i], T[i], Y_0[i],
                                           get_dYdt_epsilon, t_range,
                                           **kwargs)

        if status_i < 0:
            status = -1

        if plot_args is not None:
            if i % (len(rho)//5) == 0:
                ifig += 0.1
                plot.isotope_evolution(t_i, Y_i, particles, show=show,
                                       ifig=f"{ifig:.1f}")

        Y[i] = Y_i[:, -1]
        t[i] = t_i[-1]

    if not numpy.all(numpy.isclose(t_range[-1], t)):
        status = -1

    return Y, status

def initalize_mass_fraction(particles):
    """Get inital mass fraction from iniabu

    Parameters
    ----------
    particles : list
        List of AZN tuples of the particles in the reaction network

    Returns
    -------
    list
        List of inital mass fraction from iniabu
    """
    particles = [isotopes.AZN_to_str(particle + (1,), particle_alt=False) for particle in particles]
    return [inimf.iso[particle].abu_solar for particle in particles]

def mass_frac_to_molar_abu(particles, X):
    """Convert mass fraction to molar abundance

    Parameters
    ----------
    particles : list
        List of AZN tuples of the particles in the reaction network
    X : list
        List of mass fraction for each particle

    Returns
    -------
    list
        List of Molar abundances in [mol/g]
    """
    return [X_i/A for (A, Z), X_i in zip(particles, X)]

def molar_abu_to_mass_frac(particles, Y):
    """Convert mass fraction to molar abundance

    Parameters
    ----------
    particles : list
        List of AZN tuples of the particles in the reaction network
    Y : list
        List of Molar abundances in [mol/g]

    Returns
    -------
    list
        List of mass fraction for each particle
    """
    return [Y_i*A for (A, Z), Y_i in zip(particles, Y)]
