"""Read OPAL opacity data
"""

import numpy
from scipy.interpolate import interp1d

def read_file(fname):
    """Read OPAL Rosseland mean opacity tables

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
        file = fin.read()

    # reamove header and split into lines
    file = file.split("*** Tables ***", 1)[1]
    file = file.split("\n")[1:]

    table_len = 77

    X_i = []
    logRs = []
    logTs = []
    tables = []
    for i in range(len(file) // table_len):
        section = file[i*table_len:(i + 1)*table_len]
        #remove blank lines
        section = [line.strip() for line in section]
        section = [line for line in section if line != ""]

        header = section[0].split("=")[1:]
        header = [float(value.split()[0]) for value in header]
        X_i.append([header[0], header[2]])

        # log(R) in [g/(cm^3 K^3)]
        # To match with the standar definition of log(R) use logR + 18
        logR = section[2]
        logR = logR.split()[1:]
        logR = numpy.array([float(value) for value in logR]) - 18
        section = section[3:]

        logT = []
        table = []
        for line in section:
            line = line.replace("9.999", "Nan")
            line = line.split()
            logT.append(float(line[0]))

            line = [float(value) for value in line[1:]]
            table.append(line + [float("Nan")]*(len(logR) - len(line)))
        logT = numpy.array(logT)
        table = numpy.array(table)

        logRs.append(logR)
        logTs.append(logT)
        tables.append(table)

    if not all([(logR == logRs[0]).all() for logR in logRs[1:]]):
        raise ValueError("Tables do not share a common logR axis")

    if not all([(logT == logTs[0]).all() for logT in logTs[1:]]):
        raise ValueError("Tables do not share a common logT axis")

    logR = logRs[0]
    logT = logTs[0]
    tables = numpy.array(tables)

    # remove empty rows from Opacity Project tables
    if numpy.all(numpy.isnan(tables[:, -4:])):
        logT = logT[:-4]
        tables = tables[:, :-4]

    return X_i, logR, logT, tables

def fill_tables(opacity_data):
    """Fill missing entries from table using interpolation

    Parameters
    ----------
    opacity_data : tuple
        data from read opacity

    Returns
    -------
    tuple
        data from read opacity with out missing entries
    """
    X_i, logR, logT, tables = opacity_data
    kind = "linear"

    for i, table in enumerate(tables):
        notNan = numpy.logical_not(numpy.isnan(table))

        if numpy.any(numpy.isnan(table)):
            interp_row = [interp1d(logR[notNan[j]], table[j][notNan[j]],
                                   kind=kind, fill_value="extrapolate")(logR)
                          for j in range(len(logT))]

            interp_column = [interp1d(logT[notNan[:, j]], table[:, j][notNan[:, j]],
                                      kind=kind, fill_value="extrapolate")(logT)
                             for j in range(len(logR))]

            interp_row = numpy.array(interp_row)
            interp_column = numpy.array(interp_column).T

            tables[i] = (0.5*(interp_column + interp_row)).round(3)
    return X_i, logR, logT, tables

def regularize_tables(opacity_data):
    """Regularize opacity data to a 4-D grid

    Parameters
    ----------
    opacity_data : tuple
        Data from read_file

    Returns
    -------
    X : list
        list of hydrogen mass fractions
    Z : list
        list of metalicity
    logR : list
        axis values in the form log10(R) in [g/(cm^3 K^3)], where
        R = rho/T^3. The standard definition of logR uses R = rho/T6^3
        to convert logR into the standard version use logR + 18
    logT : list
        axis value in the form log10(T) in [k]
    tables : list
        list of opacity tables as log10(k_R) in [cm^2/g]
    """
    X_i, logR, logT, tables = fill_tables(opacity_data)

    # find metalicity values
    Z = []
    x_0 = X_i[0][0]
    for x, z in X_i:
        if x == x_0:
            Z.append(z)
        else:
            break

    # find hydrogen mass fractions
    X = []
    for i in range(0, len(X_i), len(Z)):
        if X_i[i][0] == X_i[i + len(Z) - 1][0]:
            X.append(X_i[i][0])
        else:
            break

    # check the percentage of missing tables
    data_fraction = 0
    for x in X:
        for z in Z:
            if [x, z] in X_i:
                data_fraction += 1
    data_fraction = data_fraction / (len(X)*len(Z))
    if (data_fraction < 0.90):
        print(f"WARNING: {1 - data_fraction:.2f}% of opacity tables are missing.")

    # rearange table and interpolate
    regularized_data = []
    for x in X:
        const_X = []
        for z in Z:
            if [x, z] in X_i:
                i = X_i.index([x, z])
                table = tables[i]
            else: # fill missing tables
                #TODO interpolate missing table
                raise NotImplementedError("TODO interpolate missing table")
            const_X.append(table)
        regularized_data.append(const_X)
    regularized_data = numpy.array(regularized_data)
    return X, Z, logT, logR, regularized_data
