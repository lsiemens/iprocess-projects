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
import requests
import time

from . import isotopes

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

erg_per_Mev = 1.60217663e-6 # in [erg/MeV]

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
        The energy released by the reaction in [erg]
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

        Q_value = float(head[-1])*erg_per_Mev # convert MeV to erg

        width = 13
        values  = [file[i*3 + 1][width*j:width*(j + 1)] for j in range(4)]
        values += [file[i*3 + 2][width*j:width*(j + 1)] for j in range(3)]
        values  = [float(value) for value in values]
        a.append(values)
    a = numpy.array(a).T

    def reaction_rate(T):
        """Nuclear reaction rate

        Parameters
        ----------
        T : float, array
            Tempurature in [K]

        Returns
        -------
        float, array
            Reaction rate, in [1/s] for unary rates, in [cm^3/(s mol)]
            for binary rates, in [cm^6/(s mol^2)] for trinary rates ...
        """
        T9 = 1e-9*T # Tempurature in [GK]
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

def rate_factor(Y, N, reaction_rate, rho):
    """Get rate factor

    Y_A^na Y_B^nb rho^(na + nb - 1) lambda / na! nb!

    Parameters
    ----------
    Y : array
        List of molar abundances in [mol/g]
    N : array
        number of particles involved for each type
    reaction_rate : float
        The reaction rate, in [1/s] for unary rates, in [cm^3/(s mol)]
        for binary rates, in [cm^6/(s mol^2)] for trinary rates ...
    rho : float
        density in [g/cm^3]

    Returns
    -------
    float
        Rate factor in [mol/(s g)]
    """
    Y = numpy.asarray(Y)
    N = numpy.asarray(N)
    return numpy.product(Y**N)*rho**(numpy.sum(N) - 1)*reaction_rate/numpy.product(scipy.special.factorial(N))

def make_reaction_list(dir="./nuclear/reactions/"):
    """Make reaction list

    helper file used to match reaction names with file names

    Parameters
    ----------
    dir : str, optional
        Directory to scan for nuclear reaction rates. The default is
        "./nuclear/reactions/"
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
                    path = path.replace(dir, "", 1)
                    if path.startswith("/"):
                        path = path[1:]
                    text.append(f"{reaction} {path}")
                except:
                    print(f"Failed to read \"{path}\", skipping")
                    continue

    with open(os.path.join(dir, fname_reactions), "w") as fout:
        fout.write("\n".join(text))

def read_reaction_list(dir="./nuclear/reactions/"):
    """Read reaction list

    Parameters
    ----------
    dir : str, optional
        Directory with reaction_list file. The default is
        "./nuclear/reactions/"

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

def download_reaction(reaction, file, use_set=True):
    """Download nuclear reaction rate data from JINA Reaclib

    Parameters
    ----------
    reaction : string
        Reaction equation
    file : string
        File name for the reaction rate data. It must have the format
        "path/A-BC-D-E" where A is the primary reactant, B and C are
        the seconday reactants and products respectivly, D is the
        primary product and E is the dataset label
    use_set : bool, optional
        If true, specify dataset label when downloading data from JINA
        Reaclib. The default is True
    """
    timeout = 2.5 # wait before making request from JINA Reaclib
    path, file = os.path.split(file)
    reaction = isotopes.str_to_reaction(reaction)
    num_reactants = sum([N for A, Z, N in reaction["reactants"]])
    num_products = sum([N for A, Z, N in reaction["products"]])
    format = (num_reactants, num_products)

    # split file name into A, BC, D and E. Determine how to split BC to
    # get B and C, and recombine the parts into the format expected by
    # JINA Reaclib "A(B,C)D/E"
    head, middle, tail, set = file.lower().split("-", 3)
    if format == (1, 1):
        middle = ","
    elif (format[0] == 1):
        if middle.startswith("g"):
            middle = middle[0] + "," + middle[1:]
        else:
            middle = "," + middle
    elif (format[1] == 1):
        if middle.endswith("g"):
            middle = middle[:-1] + "," + middle[-1]
        else:
            middle = middle + ","
    else:
        if middle.isalpha():
            #all one character particle labels
            middle = middle[:format[0] - 1] + "," + middle[format[0] - 1:]
        else:
            species = []
            string = ""
            for char in middle[::-1]:
                string = char + string
                if string.isalpha():
                    species = [string] + species
                    string = ""

                symbol = "".join([c for c in string if not c.isdigit()])
                if symbol in isotopes.elements:
                    species = [string] + species
                    string = ""
            species = species[:format[0] - 1] + [","] + species[format[0] - 1:]
            middle = "".join(species)
    reaction_name = f"{head}({middle}){tail}/{set}"
    if not use_set:
        reaction_name = f"{head}({middle}){tail}"

    # get rateID from JINA Reaclib
    time.sleep(timeout)
    JINA_Reaclib = "https://reaclib.jinaweb.org/"
    url_metadata = os.path.join(JINA_Reaclib, reaction_name)
    response = requests.post(url_metadata)
    response.raise_for_status()
    response = str(response.content).split("&")
    response = [line for line in response if "rateID" in line]
    rateID = response[0].split("=")[-1]

    # get reaclib file from JINA Reaclib
    time.sleep(timeout)
    php_request = f"difout.php?action=cfreaclib&rateID={rateID}&filename={file}&no910=0"
    url_data = os.path.join(JINA_Reaclib, php_request)
    response = requests.post(url_data)
    response.raise_for_status()
    response = response.content

    os.makedirs(path, exist_ok=True)
    with open(os.path.join(path, file), "wb") as fout:
        fout.write(response)
