"""Code to read JINA Reaclib reaction rates

Read nuclear reaction rates from JINA reaclib 1/2 files. The reaction
rates are stored as sets of parameters for the seven parameter function.
The resulting rate is used in a set of differential equation for the
evolution of the reactant aboundance. Details can be found at [1]_

Note the units of the reaction rate depend on the type of reaction. For
one particle reactions the units are [s^-1]. For two particle reactions
the units are [cm^3 s^-1 mol^-1] ...

Refrences
---------
[1] Cyburt, Richard H., et al. "The JINA REACLIB database: its recent
    updates and impact on type-I X-ray bursts." The Astrophysical
    Journal Supplement Series 189.1 (2010): 240.
"""


import os
import numpy
import scipy.special

import isotopes

# Reaclib chapters in the format {chapter:(# reactant, # products)}
chapters = { 1:(1, 1),
             2:(1, 2),
             3:(1, 3),
             4:(2, 1),
             5:(2, 2),
             6:(2, 3),
             7:(2, 4),
             8:(3, 1),
             9:(3, 2),
            10:(4, 2),
            11:(1, 4)}

fname_reactions = "reaction_list"

def read_file(fname):
    """Read files using REACLIB 1/2 format

    The REACLIB 1/2 file format is documented at the address
    "https://reaclib.jinaweb.org/help.php".

    Parameters
    ----------
    fname : string
        path to the file

    Returns
    -------
    reaction : dict
        Dictonary of reactants and products
    Q_value : float
        The energy released by the reaction in Mev
    reaction_rate : function
        The reaction rate function
    """
    with open(fname, "r") as fin:
        file = fin.read()
    file = file.split("\n")
    file = [line for line in file if line.strip() != ""]

    # read header
    chapter, file = int(file[0]), file[1:]
    num_input, num_output = chapters[chapter]
    num = num_input + num_output

    # remove extra header lines in REACLIB 2 files
    file = [line for line in file if len(line.strip()) > 1]

    # read sets
    a = []
    for i in range(len(file)//3):
        head = file[i*3].split()

        sort_key = lambda AZN:AZN[0]

        # get names in input/output, remove duplicates, count duplicates
        # and append number to front, convert to AZN and sort by mass number
        reactant_names = [str for str in head[:num_input]]
        reactants = list(set(reactant_names))
        reactants = [str(reactant_names.count(name)) + name for name in reactants]
        reactants = [isotopes.str_to_AZN(name) for name in reactants]
        reactants = sorted(reactants, key=sort_key, reverse=True)

        product_names = [str for str in head[num_input:num]]
        products = list(set(product_names))
        products = [str(product_names.count(name)) + name for name in products]
        products = [isotopes.str_to_AZN(name) for name in products]
        products = sorted(products, key=sort_key, reverse=True)

        reaction = {"reactants":reactants,
                    "products":products}

        Q_value = float(head[-1])

        width = 13
        values  = [file[i*3 + 1][width*j:width*(j + 1)] for j in range(4)]
        values += [file[i*3 + 2][width*j:width*(j + 1)] for j in range(3)]
        values  = [float(value) for value in values]
        a.append(values)
    a = numpy.array(a).T

    def reaction_rate(T9):
        """Nuclear reaction rate

        Parameters
        ----------
        T9 : float, array
            Tempurature in Gigakelvin

        Returns
        -------
        float, array
            Reaction rate
        """
        try:
            len(T9)
            return numpy.sum(numpy.exp(
                          a[0, :, numpy.newaxis]
                        + a[1, :, numpy.newaxis]/T9
                        + a[2, :, numpy.newaxis]/T9**(1/3)
                        + a[3, :, numpy.newaxis]*T9**(1/3)
                        + a[4, :, numpy.newaxis]*T9
                        + a[5, :, numpy.newaxis]*T9**(5/3)
                        + a[6, :, numpy.newaxis]*numpy.log(T9)), axis=0)
        except TypeError:
            return numpy.sum(numpy.exp(a[0]
                                     + a[1]/T9
                                     + a[2]/T9**(1/3)
                                     + a[3]*T9**(1/3)
                                     + a[4]*T9
                                     + a[5]*T9**(5/3)
                                     + a[6]*numpy.log(T9)), axis=0)
    return reaction, Q_value, reaction_rate

def rate_factor(Y, N, reaction_rate, T9, rho): #TODO
    """Get rate factor

    Y_A^na Y_B^nb rho^(na + nb - 1) lambda / na! nb!

    Parameters
    ----------
    reaction : dict
        Dictonary of reactants and products
    reaction_rate : function
        The reaction rate function
    T9 : float, array
            Tempurature in Gigakelvin

    Returns
    -------
    float
        Rate factor
    """
    Y = numpy.asarray(Y)
    N = numpy.asarray(N)
    return numpy.product(Y**N)*rho**(numpy.sum(N) - 1)*reaction_rate(T9)/numpy.product(scipy.special.factorial(N))

def make_reaction_list(dir="./reactions/"):
    """Make reaction list

    helper file used to match reaction names with file names

    Parameters
    ----------
    dir : str, optional
        Directory to scan for nuclear reaction rates. The default is
        "./reactions/"
    """

    text = ["# reaction file"]
    for root, dirs, files in os.walk(dir):
        for file in files:
            path = os.path.join(root, file)
            if os.path.isfile(path):
                if file.endswith(".md"):
                    continue
                if file == fname_reactions:
                    continue

                try:
                    reaction, _, _ = read_file(path)
                    reaction = isotopes.reaction_to_str(reaction)
                    text.append(f"{reaction} {path.replace(dir, '', 1)}")
                except:
                    print(f"Failed to read \"{path}\", skipping")
                    continue

    with open(os.path.join(dir, fname_reactions), "w") as fout:
        fout.write("\n".join(text))

def read_reaction_list(dir="./reactions/"):
    """Read reaction list

    Parameters
    ----------
    dir : str, optional
        Directory with reaction_list file. The default is "./reactions/"

    Returns
    -------
    dict
        Dictonary where the keys are reaction names and the value is the
        path to a reaction rate file for said reaction.
    """
    try:
        with open(os.path.join(dir, fname_reactions), "r") as fin:
            text = fin.read().split("\n")[1:]
    except FileNotFoundError:
        print("ERROR: reaction list file not found.")
        regenerate = input("Regenerate reaction list? [y/n]: ").lower()
        if regenerate.startswith("y"):
            make_reaction_list(dir)
            return read_reaction_list(dir)
        raise

    reactions_path = {}
    for line in text:
        reaction, file = line.rsplit(" ", 1)
        reactions_path[reaction] = os.path.join(dir, file)

    return reactions_path
