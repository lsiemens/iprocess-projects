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
    X_i : list
        A list of X, Z mass fraction for each table
    logR : list
        axis values in the form log10(R) in [g/(cm^3 K^3)], where
        R = rho/T^3. The standard definition of logR uses R = rho/T6^3
        to convert logR into the standard version use logR + 18
    logT : list
        axis value in the form log10(T) in [k]
    tables : list
        list of opacity tables as log10(k_R) in [cm^2/g]
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
#    print(density)

    tables = []
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
            T.append(float(line[0])*1e6)
            P.append(float(line[2])*dyne_cm2_per_MBar)
            gamma.append(float(line[8]))
        T = numpy.array(T)
        P = numpy.array(P)
        gamma = numpy.array(gamma)

        tables.append((T, P, gamma))

    return density, tables

def fill_tables(eos_data):
    density, tables = eos_data
    T = tables[0][0]

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
    regularized_data = fill_tables(eos_data)
    logRho_independent, logT, table_logP, table_gamma_rho = regularized_data

    logP = sorted(table_logP[:, 0].tolist() + table_logP[:, -1].tolist())
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

#    print(table_logRho.shape, logP.shape, logT.shape)

    return logP, logT, table_logRho, table_gamma

"""density, tables = read_file("./tables/IEOS02z8x")

lnrho, lnT, lnPs, gammas = fill_tables((density, tables))

lP, lT, lRho, lgamma = rearange((lnrho, lnT, lnPs, gammas))

print(lP.shape, lT.shape, lRho.shape, lgamma.shape)"""

#pyplot.imshow(lnPs[:, ::-1])
#pyplot.show()

"""for P_rho in lnPs:
    pyplot.plot(10**lnT, 10**P_rho)
pyplot.show()

for P_T in lnPs.T:
    pyplot.plot(10**lnrho, 10**P_T)
pyplot.show()

for gamma_rho in gammas:
    pyplot.plot(10**lnT, gamma_rho)
pyplot.show()

for gamma_T in gammas.T:
    pyplot.plot(10**lnrho, gamma_T)
pyplot.show()

print("rearange")

for Rho_P in lRho:
    pyplot.plot(10**lP, 10**Rho_P)
pyplot.show()

for Rho_T in lRho.T:
    pyplot.plot(10**lT, 10**Rho_T)
pyplot.show()

for gamma_T in lgamma.T:
    pyplot.plot(10**lT, gamma_T)
pyplot.show()

for gamma_P in lgamma:
    pyplot.plot(10**lP, gamma_P)
pyplot.show()
"""
