"""Read OPAL opacity data
"""

import numpy

def read_file(fname):
    """Read OPAL Rosseland mean opacity tables

    Parameters
    ----------
    fname : string
        file name

    Returns
    -------
    X_i : list
        A list of X, Y, Z mass fraction for each table
    logR : list
        axis values in the form log10(R), where R = rho/T6^3
    logT : list
        axis value in the form log10(T) [k]
    tables : list
        list of opacity tables as log10(k_R) [cm^2/g]
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
        X_i.append(header[:3])

        logR = section[2]
        logR = logR.split()[1:]
        logR = numpy.array([float(value) for value in logR])
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

    return X_i, logR, logT, tables

def linear_opacity(opacity_data, table_num=72):
    """Get opacity function using linear interpolation

    Parameters
    ----------
    opacity_data : tuple
        Data from read_file
    table_num : integer, optional
        Index of the table to use, this defines the composition. The
        default is 72

    Returns
    -------
    function
        Interpolated opacity function
    """
    X_i, logR, logT, tables = opacity_data
    X_i, table = X_i[table_num], tables[table_num]

    def opacity(rho, T, X, Z):
        """Interpolate opacity along rho, T, X, Z radial profiles

        Parameters
        ----------
        rho : array
            Array of densities in [g/cm^3]
        T : array
            Array of tempuratures in [K]
        X : array
            Array of hydrogen mass fraction
        Z : array
            Array of metalisity

        Returns
        -------
        array
            Opacity in [cm^2/g]
        """
        rho = numpy.asarray(rho) # density in [g/cm^3]
        T = numpy.asarray(T) # tempurature in [K]
        T6 = T/1e6 # tempurature in [MK]
        LogT = numpy.log10(T)
        LogR = numpy.log10(rho/T6**3)

        interp_fixed_T = numpy.array([numpy.interp(LogR, logR, fixed_T) for fixed_T in table])
        interp_T = numpy.array([numpy.interp(logt, logT, fixed_R) for fixed_R, logt in zip(interp_fixed_T.T, LogT)])
        return 10**interp_T # opacity in [cm^2/g]
    return opacity
