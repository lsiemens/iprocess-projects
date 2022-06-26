"""Compute opacity
"""

import numpy
import scipy.interpolate

from . import opal

def opal_opacity(fname):
    """Read opacity from OPAL Rosseland mean opacity tables

    Parameters
    ----------
    fname : string
        file name

    Returns
    -------
    function
        Opacity function
    """
    data = opal.regularize_tables(opal.read_file(fname))
    return get_opacity(data)

def get_opacity(regularized):
    """Get opacity function from regularized data

    Parameters
    ----------
    regularized : tuple
        Tuple of regularized opacity data

    Returns
    -------
    function
        Opacity function
    """
    table_X, table_Z, table_logT, table_logR, regularized_data = regularized
    def opacity(rho, T, X, Z):
        """Opacity from OPAL data table

        All input parameters must have the same shape

        Parameters
        ----------
        rho : float or array
            Density in [g/cm^3]
        T : float or array
            Temperature in [K]
        X : float or array
            Hydrogen mass fraction
        Z : float or array
            Metalisity

        Returns
        -------
        array
            Rosseland mean opacity in [cm^2/g]
        """
        logT = numpy.log10(T)
        logR = numpy.log10(rho/T**3)
        table_grid = (table_X, table_Z, table_logT, table_logR)
        parameters = numpy.array([X, Z, logT, logR]).T
        return 10**scipy.interpolate.interpn(table_grid, regularized_data, parameters, bounds_error=False, fill_value=None).T
    return opacity

def analytic_opacity(rho, T, X, Z):
    """Analytic approximation of opacity
    from https://www.astro.princeton.edu/~gk/A403/opac.pdf

    Parameters
    ----------
    rho : float or array
        Density in [g/cm^3]
    T : float or array
        Temperature in [K]
    X : float or array
        Hydrogen mass fraction
    Z : float or array
        Metalisity

    Returns
    -------
    array
        Rosseland mean opacity in [cm^2/g]
    """
    Z_star = 2*(1 + X)/(1 + 3*X - 3*Z/4)

    # electron scattering
    k_e = 0.2*(1 + X)/((1 + 2.7e11*rho/T**2)*(1 + (T/4.5e8)**0.86))
    # electron trasitions
    k_K = 4e25*(1 + X)*(Z + 0.001)*rho/T**3.5
    # H- scattering
    k_H = 1.1e-25*(Z*rho)**0.5*T**7.7
    # molecular scattering
    k_m = 0.1*Z

    # electron conductivity
    k_cond = 2.6e-7*Z_star*(T/rho)**2*(1 + (rho/2e6)**(2/3))

    # total radiative scattering
    k_rad = k_m + 1/(1/k_H +1/(k_e + k_K))

    # total effective opacity
    return 1/(1/k_rad + 1/k_cond)
