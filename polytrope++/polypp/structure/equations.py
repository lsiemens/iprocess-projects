"""Solve the equation of stellar structure
"""

import numpy
import scipy.integrate

from polypp.constants import G, R, sigma, M_sol, R_sol, L_sol

def get_dPrLTdM(get_dYdt_epsilon, X_metalicity, opacity, EOS, Y_interpolate):
    """Get dPrLTdM structure equations

    Parameters
    ----------
    get_dYdt_epsilon : function
        Function that gets dYdt and epsilon at defined density and
        tempurature for a given nuclear network
    X_metalicity : function
        Function that returns the total hydrogen mass fraction and
        the metalicity
    opacity : function
        The opacity interpolation function
    EOS : function
        The equation of state function
    Y_interpolate : function
        Function that interpolates the molar abundance vector Y
        as a function of mass coorinate

    Returns
    -------
    function
        Function dPrLTdM(M, PrLT), is the mass coordinate derivative of
        the structure variables PrLT, where the variables are Pressure
        in [dyne/cm^2], radius in [cm], Luminosity in [erg/s] and
        Tempurature in [K]
    """
    def dPrLTdM(M, PrLT):
        """stellar structure differental equation

        dPrLTdM(M, PrLT), is the mass coordinate derivative of the
        structure variables PrLT, where the variables are Pressure in
        [dyne/cm^2], radius in [cm], Luminosity in [erg/s] and
        Tempurature in [K]

        Parameters
        ----------
        M : float
            Mass coordinate in [g]
        PrLT : list
            the structure variables PrLT, where the variables are Pressure
            in [dyne/cm^2], radius in [cm], Luminosity in [erg/s] and
            Tempurature in [K]

        Returns
        -------
        list
            The differential dPrLTdM in [dyne/(cm^2 g), cm/g, erg/(s g), K/g]
        """
        P, r, L, T = numpy.abs(PrLT)

        Y = Y_interpolate(M) # molar abundance vector in [mol/g]
        X, Z = X_metalicity(Y) # hydrogen mass fraction and metalicity

        rho, gamma = EOS(P, T, X, Z) # with rho in [g/cm^3] and gamma is unitless

        # nuclear energy production density epsilon
        _, epsilon_fun = get_dYdt_epsilon(rho, T)
        # use Y molar abundance
        epsilon = epsilon_fun(0, Y)

        # rosseland mean opacity
        k_R = opacity(rho, T, X, Z)[0]

        # structure equations
        dPdM = -G*M/(4*numpy.pi*r**4) # hydrostatic equilbrium
        drdM = 1/(4*numpy.pi*r**2*rho) # radial distribution
        dLdM = epsilon # energy production
        dTdM_rad = -3*k_R*L/(256*numpy.pi**2*r**4*sigma*T**3) # thermal equilibrium

        # Schwarzschild stability criterion
        # TODO consider rearanging terms to remove division by zero
        if (P/T)*dTdM_rad/dPdM < (1 - 1/gamma):
            dTdM = dTdM_rad
        else:
            dTdM = (1 - 1/gamma)*(T/P)*dPdM

        return [dPdM, drdM, dLdM, dTdM]

    def radiative_cells(M, PrLT):
        """Determine if a cell is radiative

        Parameters
        ----------
        M : float
            Mass coordinate in [g]
        PrLT : list
            the structure variables PrLT, where the variables are Pressure
            in [dyne/cm^2], radius in [cm], Luminosity in [erg/s] and
            Tempurature in [K]

        Returns
        -------
        array
            boolean array that is true where the cell is radiative
        """
        P, r, L, T = numpy.abs(PrLT)

        Y = Y_interpolate(M) # molar abundance vector in [mol/g]
        X, Z = X_metalicity(Y) # hydrogen mass fraction and metalicity

        rho, gamma = EOS(P, T, X, Z) # with rho in [g/cm^3] and gamma is unitless

        # rosseland mean opacity
        k_R = opacity(rho, T, X, Z)

        # structure equations
        dPdM = -G*M/(4*numpy.pi*r**4) # hydrostatic equilbrium
        dTdM_rad = -3*k_R*L/(256*numpy.pi**2*r**4*sigma*T**3) # thermal equilibrium

        # Schwarzschild stability criterion
        # TODO consider rearanging terms to remove division by zero
        return (P/T)*dTdM_rad/dPdM < (1 - 1/gamma)

    def convective_zones(M, PrLT):
        """Get list of convective zones

        Parameters
        ----------
        M : float
            Mass coordinate in [g]
        PrLT : list
            the structure variables PrLT, where the variables are Pressure
            in [dyne/cm^2], radius in [cm], Luminosity in [erg/s] and
            Tempurature in [K]

        Returns
        -------
        list
            list of boolean arrays that are true where the cell is convective.
            Each element of the list contains one continous convective zone
        """
        convective = numpy.logical_not(radiative_cells(M, PrLT))
        bands = []

        isBand = False # currently in convective band
        for i, zone in enumerate(convective):
            if (not isBand) and zone:
                isBand = True
                bands.append([False]*len(convective))

            if isBand:
                # add point to current active convective zone
                # and stop at the next non-convective zone
                bands[-1][i] = zone
                isBand = zone
        return numpy.array(bands)

    def convective_mixing(M, PrLT):
        """Get list of convective zones

        Parameters
        ----------
        M : float
            Mass coordinate in [g]
        PrLT : list
            the structure variables PrLT, where the variables are Pressure
            in [dyne/cm^2], radius in [cm], Luminosity in [erg/s] and
            Tempurature in [K]

        Returns
        -------
        array
            The molar abundance vector in [mol/g] after convective mixing
        """
        bands = convective_zones(M, PrLT)
        Y = Y_interpolate(M) # molar abundance vector in [mol/g]

        delta_M = numpy.diff(M.tolist() + [0])

        # replace Y values in each band with the weighted mean
        for band in bands:
            Y_mixed = numpy.sum(Y[:, band]*delta_M[band]/numpy.sum(delta_M[band]), axis=1)
            Y[:, band] = Y_mixed[:, numpy.newaxis]
        return Y
    return dPrLTdM, convective_mixing

def zP(M, PrLT):
    """Zero pressure condition, terminate solution at zero pressure

    Parameters
    ----------
    M : float
        Mass coordinate in [g]
    PrLT : list
        the structure variables PrLT, where the variables are Pressure
        in [dyne/cm^2], radius in [cm], Luminosity in [erg/s] and
        Tempurature in [K]

    Returns
    -------
    float
        Pressure
    """
    P, r, L, T = PrLT
    return P
zP.terminal = True

def zr(M, PrLT):
    """Zero radius condition, terminate solution at zero radius

    Parameters
    ----------
    M : float
        Mass coordinate in [g]
    PrLT : list
        the structure variables PrLT, where the variables are Pressure
        in [dyne/cm^2], radius in [cm], Luminosity in [erg/s] and
        Tempurature in [K]

    Returns
    -------
    float
        Radius
    """
    P, r, L, T = PrLT
    return r
zr.terminal = True

def zL(M, PrLT):
    """Zero luminosity condition, terminate solution at zero luminosity

    Parameters
    ----------
    M : float
        Mass coordinate in [g]
    PrLT : list
        the structure variables PrLT, where the variables are Pressure
        in [dyne/cm^2], radius in [cm], Luminosity in [erg/s] and
        Tempurature in [K]

    Returns
    -------
    float
        Luminosity
    """
    P, r, L, T = PrLT
    return L
zL.terminal = True

def zT(M, PrLT):
    """Zero tempurature condition, terminate solution at zero tempuratue

    Parameters
    ----------
    M : float
        Mass coordinate in [g]
    PrLT : list
        the structure variables PrLT, where the variables are Pressure
        in [dyne/cm^2], radius in [cm], Luminosity in [erg/s] and
        Tempurature in [K]

    Returns
    -------
    float
        Tempurature
    """
    P, r, L, T = PrLT
    return T
zT.terminal = True

def core_error(M, PrLT):
    """Boundary condition error at the core

    Parameters
    ----------
    M : float
        Mass coordinate in [g]
    PrLT : list
        the structure variables PrLT, where the variables are Pressure
        in [dyne/cm^2], radius in [cm], Luminosity in [erg/s] and
        Tempurature in [K]

    Returns
    -------
    list
        List of error terms that should be zero at the core
    """
    P_s, r_s, L_s, T_s = PrLT[:, 0]
    P_c, r_c, L_c, T_c = PrLT[:, -1]
    M_s = M[0]
    M_c = M[-1]
    return [M_c/M_s, r_c/r_s, L_c/L_s]

def solve_core(dPrLTdM, M, r, T, **kwargs):
    """Solve stellar structure equations given surface parameters

    Parameters
    ----------
    dPrLTdM : function
        function that computes the differential dPrLTdM in
        [dyne/(cm^2 g), cm/g, erg/(s g), K/g]
    M : float
        Mass coordinate in [g]
    r : float
        radius in [cm]
    T : float
        surface tempurature in [K]
    kwargs : dictonary
        arguments for scipy.integrate.solve_ivp

    Returns
    -------
    M : float
        Mass coordinate in [g]
    PrLT : list
        the structure variables PrLT, where the variables are Pressure
        in [dyne/cm^2], radius in [cm], Luminosity in [erg/s] and
        Tempurature in [K]
    status : bool
        Status of solve_ivp solution
    error : float
        The norm of the solutions core_error
    """
    k_R = 0.2*(1 + 0.7)

    P = (2/3)*G*M/(r**2*k_R)
    L = 4*numpy.pi*r**2*sigma*T**4

    events = [zr, zL]
    PrLT = [P, r, L, T]
    solve = scipy.integrate.solve_ivp(dPrLTdM, (M, 0), PrLT, events=events, method="Radau", **kwargs)
    if not solve.success:
        print(solve.message)

    M, PrLT = solve["t"], solve["y"]
    error = core_error(M, PrLT)
    return M, PrLT, solve.status, numpy.linalg.norm(error)

def estimate_rT(M):
    """Estimate radius and tempurature

    Parameters
    ----------
    M : float
        Mass in [g]

    Returns
    -------
    r : float
        radius in [cm]
    T : float
        surface tempurature in [K]
    """
    R = R_sol*(M/M_sol)**0.8
    L = L_sol*(M/M_sol)**3.5
    T = (L/(4*numpy.pi*R**2*sigma))**(1/4)
    return R, T

def minimize_grid(dPrLTdM, M, r_range, T_range, n=5, return_values=False):
    """Solve stellar structure equations given surface parameters

    Parameters
    ----------
    dPrLTdM : function
        function that computes the differential dPrLTdM in
        [dyne/(cm^2 g), cm/g, erg/(s g), K/g]
    M : float
        Mass coordinate in [g]
    r_range : tuple
        radius range in [cm]
    T_range : tuple
        surface tempurature range in [K]
    n : integer, optional
        resolution of parameter grid. The default is 5
    return_values : array
        If True then return the array of error values

    Returns
    -------
    min_param : tuple
        Tuple of the radius, tempurature pair which minimized the core error
    error : float
        The minimum error found on the parameter grid
    values : array
        The array of computed core errors. This parameter is only returned
        if return_values is True
    """
    min_error = float("inf")
    min_param = None
    k = 0
    if return_values:
        values = numpy.zeros(shape=(n,n))
    for i, r in enumerate(numpy.linspace(*r_range, n)):
        for j, T in enumerate(numpy.linspace(*T_range, n)):
            k += 1
            if k % max(1, n**2//10) == 0:
                print(k)

            _, _, _, error = solve_core(dPrLTdM, M, r, T)
            values[i, j] = error
            if error < min_error:
                min_error = error
                min_param = (r, T)

    if return_values:
        return min_param, min_error, values
    else:
        return min_param, min_error
