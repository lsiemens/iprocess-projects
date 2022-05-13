"""Solve Lane-Emden equation transformed to use mass coordinates
"""

import numpy
import scipy.integrate

G = 6.674e-8 # in [dyne cm^2 / g^2]
M_sol = 1.988e33 # in [g]
R_sol = 6.957e10 # in [cm]


# define y = <theta, chi'>, t = phi
def get_dydt(alpha, beta, n):
    """Get dydt(t, y)

    let t = phi and y = <theta, chi>

    Parameter
    ---------
    n : float
        Polytropic index

    Returns
    -------
    f(t, y)
        The function dydt
    """
    def dydt(t, y):
        """a system of first order equations

        Parameters
        ----------
        t : float
            Equivelent to the variable M
        y : array
            The vector <theta, chi>
        """
        theta, chi = y
        dtheta_dt = -alpha*t/((n + 1)*chi**4*theta**n)
        dchi_dt = beta/(chi**2*theta**n)
        return [dtheta_dt, dchi_dt]
    return dydt

def solve_polytrope(core_pressure, core_density, n):
    """Find polytrope

    Parameters
    ----------
    core_pressure : float
        Pressure at the stars core in
    core_density : float
        Density at the stars core in [g/cm^3]
    n : float
        The polytropic index n, n = 1 convective stars, n = 3 radiative
        zone of main sequence star.
    """

    alpha = G*M_sol**2/(4*numpy.pi*R_sol**4*core_pressure)
    beta = M_sol/(4*numpy.pi*R_sol**3*core_density)

    t_range = (0, 100)
    surface_event = lambda phi, y:y[0] # stop simulation when theta = 0
    surface_event.terminal = True

    #[core_density, radius ~ 0]
    y0 = [1.0, 1e-3]

    solution = scipy.integrate.solve_ivp(get_dydt(alpha, beta, n), t_range, y0, events=surface_event)
    phi = solution["t"]
    theta, chi = solution["y"]

    return M_sol*phi, core_density*theta**n, R_sol*chi
