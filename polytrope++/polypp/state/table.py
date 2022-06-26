import numpy
import scipy.interpolate

from . import opal_eos
from polypp.constants import R

def opal_data(fname):
    """Load equation of state function from OPAL EOS tables

    Parameters
    ----------
    fname : string
        file name

    Returns
    -------
    eos : tuple
        Equation of state functions (Rho(P, T, X, Z), gamma(P, T, X, Z))
    """
    data = opal_eos.rearange(opal_eos.read_file(fname))
    return get_eos(data)

def get_eos(data):
    """get equations of state functions from EOS data tables

    Parameters
    ----------
    data : tuple
        data from opal_eos.rearange

    Returns : tuple
        Equation of state functions (Rho(P, T, X, Z), gamma(P, T, X, Z))
    """
    table_logP, table_logT, table_logRho, table_gamma = data
    def eos(P, T, X, Z):
        logP = numpy.log10(P)
        logT = numpy.log10(T)
        table_grid = (table_logP, table_logT)
        parameters = numpy.array([logP, logT]).T
        Rho = 10**scipy.interpolate.interpn(table_grid, table_logRho, parameters, bounds_error=False, fill_value=None).T
        gamma = scipy.interpolate.interpn(table_grid, table_gamma, parameters, bounds_error=False, fill_value=None).T
        return Rho, gamma
    return eos

def ideal_eos(P, T, X, Z):
    """Ideal gas equation of state

    Parameters
    ----------
    P : float
        Pressure
    T : float
        Tempurature
    X : float
        hydrogen mass fraction
    Z : float
        metalicity
    """
    Rho = mu(X, Z)*P/(R*T)
    return Rho, 5/3.0 + 0*Rho

def approximate_eos(P, T, X, Z):
    pass

def mu(X, Z):
    """mean molecular weight

    Parameters
    ----------
    X : float
        hydrogen mass fraction
    Z : float
        metalicity

    Returns
    -------
    mu : float
        mean molecular weight
    """
    return 4/(5*X + 3 - Z)
