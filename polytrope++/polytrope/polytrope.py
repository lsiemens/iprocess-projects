"""Solve Lane-Emden equation transformed to use mass coordinates
"""

import numpy
import scipy.integrate

G = 6.674e-8 # in [dyne cm^2 / g^2]
M_sol = 1.988e33 # in [g]
R_sol = 6.957e10 # in [cm]
L_sol = 3.83e33 # in [erg/s]

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

def solve_polytrope(core_density, core_pressure, n, phi_range=(0, 10),
                    num_steps=None, **kwargs):
    """Find polytrope

    Parameters
    ----------
    core_density : float
        Density at the stars core in [g/cm^3]
    core_pressure : float
        Pressure at the stars core in
    n : float
        The polytropic index n, n = 1 convective stars, n = 3 radiative
        zone of main sequence star.
    phi_range : tuple, optional
        Range of phi values to solve over. The default is (0, 10)
    num_steps : integer or None, optional
        If integer, phi_range is split into num_steps equal steps. If
        None, then leave the choice upto solve_ivp. The default is None
    kwargs : dict
        Key word arguments for scipy solve_ivp

    Returns
    -------
    phi : array
        Mass coordinate array, where M = M_sol*phi
    theta : array
        polytrope parameter, where rho = core_density*theta^n
    chi : array
        Radius array, where r = R_sol*chi
    """

    alpha = G*M_sol**2/(4*numpy.pi*R_sol**4*core_pressure)
    beta = M_sol/(4*numpy.pi*R_sol**3*core_density)

    # set inital theta and chi
    # theta_0 = 1.0 ie density at the center is core_density
    # chi_0 = 1e-9 ie radius at the center is 1e-9 R_0 ~ 1 meter
    y0 = [1.0, 1e-9]

    surface_event = lambda phi, y:y[0] # stop simulation when theta = 0
    surface_event.terminal = True

    dydt = get_dydt(alpha, beta, n)

    if num_steps is not None:
        phi = numpy.linspace(*phi_range, num_steps)

        solution = scipy.integrate.solve_ivp(dydt, phi_range,
                                             y0, t_eval=phi,
                                             events=surface_event, **kwargs)
    else:
        solution = scipy.integrate.solve_ivp(dydt, phi_range, y0,
                                             events=surface_event, **kwargs)

    print(f"Message from solve_ivp: {solution.message}")
    phi = solution["t"]
    theta, chi = solution["y"]

    return phi, theta, chi

def polytrope_to_cgs(data, core_density, core_pressure, n):
    """Calculate mass, radius, density and pressure in cgs

    Parameters
    ----------
    data : tuple
        Tuple containing phi, theta and chi arrays
    core_density : float
        Core density in [g/cm^3]
    core_pressure : float
        Core pressure in [dyn/cm^2]
    n : float
        The polytropic index n

    Returns
    -------
    m : array
        Mass coordinate in [g]
    r : array
        Radius in [cm]
    rho : array
        Density in [g/cm^3]
    P : array
        Pressure in [dyn/cm^2]
    """
    phi, theta, chi = data

    m = M_sol*phi
    r = R_sol*chi
    rho = core_density*theta**n
    P = core_pressure*theta**(n + 1)

    return m, r, rho, P
