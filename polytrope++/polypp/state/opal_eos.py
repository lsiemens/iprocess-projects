"""Read OPAL opacity data
"""

dyne_cm2_per_MBar = 1e12 # in [dyne/(cm^2 MBar)]

import numpy
from scipy.interpolate import interp1d
from matplotlib import pyplot

def read_file(fname):
    """Read OPAL Equation of State data table

    Parameters
    ----------
    fname : string
        file name

    Returns
    -------
    density : list
        list of table densities
    tables : list
        tables of tempurature, pressure and gamma at the given densities
    """
    with open(fname, "r") as fin:
        file = fin.read().strip().split("\n")

    header, file = file[:2], file[2:]
    # TODO readin X, Z mass fractions

    file = [line.strip() for line in file if line.strip() != ""]
    headers = [line for line in file if "density" in line]
    table_len = []
    density = []
    for line in headers:
        header = line.split()
        table_len.append(int(header[1]))
        density.append(float(header[5]))

    density = numpy.array(density)

    tables = []

    if "density" not in file[0]:
        raise RuntimeError("Expected start of first table")

    # split columns names if seporated by double space
    table_legend = [item.strip() for item in file[1].split("  ")]
    # remove empty items
    table_legend = [item for item in table_legend if len(item) > 0]

    # get ids of colomns of interest
    T_id = table_legend.index("T6")
    P_id = table_legend.index("P(MB)")
    gamma_id = table_legend.index("g1")

    for length in table_len:
        if "density" not in file[0]:
            raise RuntimeError("Expected start of table")
        file = file[2:]
        table, file = file[:length], file[length:]
        T = []
        P = []
        gamma = []
        for line in table:
            line = line.split()
            T.append(float(line[T_id])*1e6)
            P.append(float(line[P_id])*dyne_cm2_per_MBar)
            gamma.append(float(line[gamma_id]))
        T = numpy.array(T)
        P = numpy.array(P)
        gamma = numpy.array(gamma)

        tables.append((T, P, gamma))
    return density, tables

def fill_tables(eos_data):
    """fill missing entries in table

    Parameters
    ----------
    eos_data : tuple
        eos data from read_file

    Returns
    -------
    logRho : array
        logarithm of density axis
    logT : array
        logarithm of tempurature axis
    logP : array
        logarithm of pressure table
    gamma : array
        gamma table
    """
    density, tables = eos_data
    T = tables[0][0][::-1]

    logP = []
    gamma = []
    for table in tables[::-1]:
        T_partial, P_partial, gamma_partial = table
        logP_interp = interp1d(numpy.log10(T_partial), numpy.log10(P_partial), fill_value="extrapolate", kind="linear")
        gamma_interp = interp1d(numpy.log10(T_partial), gamma_partial, fill_value="extrapolate", kind="linear")

        logP.append(logP_interp(numpy.log10(T)))
        gamma.append(gamma_interp(numpy.log10(T)))

    logP = numpy.array(logP)
    gamma = numpy.array(gamma)

    return numpy.log10(density[::-1]), numpy.log10(T), logP, gamma

def rearange(eos_data):
    """rearange tables of logP(logRho, logT) to the form logRho(logP, logT)

    Parameters
    ----------
    eos_data : tuple
        eos data from read_file

    Returns
    -------
    logP : array
        logarithm of pressure axis
    logT : array
        logarithm of tempurature axis
    logRho : array
        logarithm of density table
    gamma : array
        gamma table
    """
    regularized_data = fill_tables(eos_data)
    logRho_independent, logT, table_logP, table_gamma_rho = regularized_data

    logP = table_logP[:, 0].tolist() + table_logP[:, -1].tolist()
    # sort and remove duplicate points
    logP = sorted(list(set(logP)))
    logP = numpy.array(logP)

    table_logRho = []
    table_gamma = []
    for logP_dependent, gamma_rho in zip(table_logP.T, table_gamma_rho.T):
        # logP, logRho at fixed T. Convert dependent variable logP to independent variable
        logRho_interp = interp1d(logP_dependent, logRho_independent, fill_value="extrapolate", kind="linear")
        gamma_interp = interp1d(logP_dependent, gamma_rho, fill_value="extrapolate", kind="linear")


        table_logRho.append(logRho_interp(logP))
        table_gamma.append(gamma_interp(logP))

    table_logRho = numpy.array(table_logRho).T
    table_gamma = numpy.array(table_gamma).T

    return logP, logT, table_logRho, table_gamma
