"""Solve the equation of stellar structure
"""

import numpy

G = 6.674e-8 # gravitational constant in [dyne cm^2 / g^2]
R = 8.314e7 # gass constant in [erg/(g K)]
sigma = 5.67037442e-5 # Stefan-Boltzmann constant in [erg/(s cm^2 K^4)]

def get_dPrLTdM(get_dYdt_epsilon, X_metalicity, opacity, Y_interpolate):
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

        # mean molecular weight
        mu = 4/(5*X + 3 - Z)
        gamma = 5/3

        rho = mu*P/(R*T) # in [g/cm^3]

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
        if (P/T)*dTdM_rad/dPdM < (1 - 1/gamma):
            dTdM = dTdM_rad
        else:
            dTdM = (1 - 1/gamma)*(T/P)*dPdM

        return [dPdM, drdM, dLdM, dTdM]
    return dPrLTdM

def zP(M, PrLT):
    P, r, L, T = PrLT
    return P
zP.terminal = True

def zr(M, PrLT):
    P, r, L, T = PrLT
    return r
zr.terminal = True

def zL(M, PrLT):
    P, r, L, T = PrLT
    return L
zL.terminal = True

def zT(M, PrLT):
    P, r, L, T = PrLT
    return T
zT.terminal = True

def BC(M, PrLT):
    P_s, r_s, L_s, T_s = PrLT[:, 0]
    P_c, r_c, L_c, T_c = PrLT[:, -1]
    M_s = M[0]
    M_c = M[-1]
    return [M_c/M_s, r_c/r_s, L_c/L_s]
