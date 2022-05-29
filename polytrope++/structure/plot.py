"""Plots of stellar structure
"""

"""Plots and graphs
"""

from matplotlib import pyplot

def PrLT(M, PrLT, massCoord=True, show=True, ifig=None):
    """Plot isotope mass fraction vs coordinate (mass or radius)

    Parameters
    ----------
    M : float
        Mass coordinate in [g]
    PrLT : list
        the structure variables PrLT, where the variables are Pressure in [dyne/cm^2],
        radius in [cm], Luminosity in [erg/s] and Tempurature in [K]
    massCoord : bool, optional
        If True then x is treated as a mass coordinate, otherwise it is
        assumed to be a radial coodinate. The default is True
    show : bool, optional
        If True, show the network. The default is False
    ifig : integer, optional
        Figure number. The default is None
    """
    P, r, L, T = PrLT

    if ifig is not None:
        pyplot.close(ifig)
        pyplot.figure(ifig)

    if massCoord:
        pyplot.plot(M/M.max(), P/P.max(), label="pressure $y = P$")
        pyplot.plot(M/M.max(), r/r.max(), label="radius $y = R$")
        pyplot.plot(M/M.max(), L/L.max(), label="luminosity $y = L$")
        pyplot.plot(M/M.max(), T/T.max(), label="tempurature $y = T$")
        pyplot.xlabel("mass coordinate")
    else:
        pyplot.plot(r, P/P.max(), label="pressure $y = P$")
        pyplot.plot(r, M/M.max(), label="mass $y = M$")
        pyplot.plot(r, L/L.max(), label="luminosity $y = L$")
        pyplot.plot(r, T/T.max(), label="tempurature $y = T$")
        pyplot.xlabel("radius [cm]")
    pyplot.ylabel("parameter $y$ fraction of maximum")
    pyplot.legend()

    if show:
        if ifig is not None:
            pyplot.show(ifig)
        else:
            pyplot.show()

"""
ifig += 1; pyplot.close(ifig); pyplot.figure(ifig)
pyplot.plot(M/M.max(), r/r.max(), label="radius $y = R$")
pyplot.plot(M/M.max(), rho/rho.max(), label="density $y = \\rho$")
pyplot.plot(M/M.max(), P/P.max(), label="pressure $y = P$")
pyplot.plot(M/M.max(), T/T.max(), label="tempurature $y = T$")
pyplot.plot(M/M.max(), L/L.max(), label="luminosity $y = L$")
pyplot.plot(M/M.max(), k_R/k_R.max(), label="mean opacity $y = k_R$")
pyplot.plot(M/M.max(), power_density/power_density.max(), label="power density $y = \\varepsilon \\rho$")
pyplot.xlabel("mass coordinate")
pyplot.ylabel("parameter $y$ fraction of maximum")
pyplot.legend()
pyplot.show(ifig)

ifig += 1; pyplot.close(ifig); pyplot.figure(ifig)
pyplot.plot(r, M/M.max(), label="mass $y = M$")
pyplot.plot(r, rho/rho.max(), label="density $y = \\rho$")
pyplot.plot(r, P/P.max(), label="pressure $y = P$")
pyplot.plot(r, T/T.max(), label="tempurature $y = T$")
pyplot.plot(r, L/L.max(), label="luminosity $y = L$")
pyplot.plot(r, k_R/k_R.max(), label="mean opacity $y = k_R$")
pyplot.plot(r, power_density/power_density.max(), label="power density $y = \\varepsilon \\rho$")
pyplot.xlabel("radius [cm]")
pyplot.ylabel("parameter $y$ fraction of maximum")
pyplot.legend()
pyplot.show(ifig)

ifig += 1; pyplot.close(ifig); pyplot.figure(ifig)
pyplot.plot(r, rho)
pyplot.xlabel("radius [cm]")
pyplot.ylabel("rho [g/cm^3]")
pyplot.title("density")
pyplot.show(ifig)

ifig += 1; pyplot.close(ifig); pyplot.figure(ifig)
pyplot.plot(r, T)
pyplot.xlabel("radius [cm]")
pyplot.ylabel("temperature [K]")
pyplot.title("Temperature")
pyplot.show(ifig)

ifig += 1; pyplot.close(ifig); pyplot.figure(ifig)
pyplot.semilogy(r, mfp)
pyplot.xlabel("radius [cm]")
pyplot.ylabel("mean free path [cm]")
pyplot.title("Opacity: Mean free path of light")
pyplot.show(ifig)

ifig += 1; pyplot.close(ifig); pyplot.figure(ifig)
pyplot.plot(r, power_density)
pyplot.xlabel("radius [cm]")
pyplot.ylabel("power density [erg/(s cm^3)]")
pyplot.title("Nuclear power density")
pyplot.show(ifig)

ifig += 1; pyplot.close(ifig); pyplot.figure(ifig)
pyplot.plot(r, m)
pyplot.xlabel("radius [cm]")
pyplot.ylabel("mass [g]")
pyplot.title("total mass")
pyplot.show(ifig)
"""
