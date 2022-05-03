"""Code to read JINA Reaclib reaction rates
"""

import numpy
from numpy import sum, log, exp

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

def read_file(fname):
    """Read files using REACLIB 1/2 format

    The REACLIB 1/2 file format is documented at the address
    "https://reaclib.jinaweb.org/help.php".

    Parameters
    ----------
    fname : string
        path to the file
    """
    with open(fname, "r") as fin:
        file = fin.read()
    file = file.split("\n")
    file = [line for line in file if line.strip() != ""]

    # read header
    chapter, file = int(file[0]), file[1:]

    # read sets
    a = []
    for i in range(len(file)//3):
        head = file[i*3].split()
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
            return numpy.sum(numpy.exp(a[0]
                                     + a[1]/T9
                                     + a[2]/T9**(1/3)
                                     + a[3]*T9**(1/3)
                                     + a[4]*T9
                                     + a[5]*T9**(5/3)
                                     + a[6]*numpy.log(T9)), axis=0)
        except ValueError:
            return numpy.sum(numpy.exp(
                          a[0, :, numpy.newaxis]
                        + a[1, :, numpy.newaxis]/T9
                        + a[2, :, numpy.newaxis]/T9**(1/3)
                        + a[3, :, numpy.newaxis]*T9**(1/3)
                        + a[4, :, numpy.newaxis]*T9
                        + a[5, :, numpy.newaxis]*T9**(5/3)
                        + a[6, :, numpy.newaxis]*numpy.log(T9)), axis=0)
    return reaction_rate

nrr = read_file("n15-pg-o16-li10")
from matplotlib import pyplot

print(nrr(0.5))

T9 = numpy.linspace(1e-3, 10)
pyplot.plot(T9, nrr(T9))
pyplot.show()
