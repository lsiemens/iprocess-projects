import numpy
import scipy.interpolate

from . import opal_eos

def opal_data(fname):
    data = opal_eos.rearange(opal_eos.read_file(fname))
    return get_eos(data)

def get_eos(data):
    table_logP, table_logT, table_logRho, table_gamma
    def eos(P, T, X, Z):
        logP = numpy.log10(P)
        logT = numpy.log10(T)
        table_grid = (table_logP, table_logT)
        parameters = numpy.array([logP, logT]).T
        Rho = 10**scipy.interpolate.interpn(table_grid, table_logRho, parameters, bounds_error=False, fill_value=None).T
        gamma = scipy.interpolat.interpn(table_grid, table_gamma, parameters, bounds_error=False, fill_value=None).T
        return Rho, gamma
    return eos

def ideal_eos(P, T, X, Z):
    pass

def approximate_eos(P, T, X, Z):
    pass
