"""
General Numerical Solver for the 1D Time-Dependent Schrodinger Equation.

Authors:
- Jake Vanderplas <vanderplas@astro.washington.edu>
- Andre Xuereb (imaginary time propagation, normalized wavefunction
- Luke Siemens <luke@lsiemens.com>

For a theoretical description of the algorithm, please see
http://jakevdp.github.com/blog/2012/09/05/quantum-python/

For the original version of this code please see
http://github.com/jakevdp/pySchrodinger/

This version of schrodinger.py has been modified to produce results in
terms of momentum rather than wavenumber.

License: BSD style
Please feel free to use and modify this, but keep the above information.
"""

import numpy as np
from scipy import fftpack


class Schrodinger(object):
    """
    Class which implements a numerical solution of the time-dependent
    Schrodinger equation for an arbitrary potential
    """
    def __init__(self, x, psi_x0, V_x, p0=None, hbar=1, m=1, t0=0.0):
        """
        Parameters
        ----------
        x : array_like, float
            Length-N array of evenly spaced spatial coordinates
        psi_x0 : array_like, complex
            Length-N array of the initial wave function at time t0
        V_x : array_like, float
            Length-N array giving the potential at each x
        p0 : float
            The minimum value of p.  Note that, because of the workings
            of the Fast Fourier Transform, the momentum will be defined
            in the range
              p0 < p < 2*pi*hbar / dx ,
            where dx = x[1]-x[0].  If you expect nonzero momentum
            outside this range, you must modify the inputs accordingly.
            If not specified, p0 will be calculated such that the range
            is [-p0,p0]
        hbar : float
            Value of Planck's constant (default = 1)
        m : float
            Particle mass (default = 1)
        t0 : float
            Initial time (default = 0)
        """
        # Validation of array inputs
        self.x, psi_x0, self.V_x = map(np.asarray, (x, psi_x0, V_x))
        N = self.x.size
        assert self.x.shape == (N,)
        assert psi_x0.shape == (N,)
        assert self.V_x.shape == (N,)

        # Validate and set internal parameters
        assert hbar > 0
        assert m > 0
        self.hbar = hbar
        self.m = m
        self.t = t0
        self.dt_ = None
        self.N = len(x)
        self.dx = self.x[1] - self.x[0]
        self.dp = 2 * np.pi * hbar / (self.N * self.dx)

        # Set momentum scale
        if p0 == None:
            self.p0 = -0.5 * self.N * self.dp
        else:
            assert p0 < 0
            self.p0 = p0
        self.p = self.p0 + self.dp * np.arange(self.N)

        self.psi_x = psi_x0
        self.compute_p_from_x()

        # Variables which hold steps in evolution
        self.x_evolve_half = None
        self.x_evolve = None
        self.p_evolve = None

    def _set_psi_x(self, psi_x):
        assert psi_x.shape == self.x.shape
        self.psi_mod_x = (psi_x * np.exp(-1j * self.p[0] * self.x / self.hbar)
                          * self.dx / np.sqrt(2 * np.pi * self.hbar))
        self.psi_mod_x /= self.norm
        self.compute_p_from_x()

    def _get_psi_x(self):
        return (self.psi_mod_x * np.exp(1j * self.p[0] * self.x / self.hbar)
                * np.sqrt(2 * np.pi * self.hbar) / self.dx)

    def _set_psi_p(self, psi_p):
        assert psi_p.shape == self.x.shape
        self.psi_mod_p = psi_p * np.exp(1j * self.x[0] * self.dp
                                        * np.arange(self.N) / self.hbar)
        self.compute_x_from_p()

    def _get_psi_p(self):
        return self.psi_mod_p * np.exp(-1j * self.x[0] * self.dp
                                        * np.arange(self.N) / self.hbar)

    def _get_dt(self):
        return self.dt_

    def _set_dt(self, dt):
        assert dt != 0
        if dt != self.dt_:
            self.dt_ = dt
            self.x_evolve_half = np.exp(-0.5 * 1j * self.V_x * self.dt
                                         / self.hbar)
            self.x_evolve = self.x_evolve_half ** 2
            self.p_evolve = np.exp(-0.5 * 1j * (self.p ** 2) * self.dt
                                    / (self.m * self.hbar))

    def _get_norm(self):
        return self.wf_norm(self.psi_mod_x)

    psi_x = property(_get_psi_x, _set_psi_x)
    psi_p = property(_get_psi_p, _set_psi_p)
    norm = property(_get_norm)
    dt = property(_get_dt, _set_dt)

    def compute_p_from_x(self):
        self.psi_mod_p = fftpack.fft(self.psi_mod_x)

    def compute_x_from_p(self):
        self.psi_mod_x = fftpack.ifft(self.psi_mod_p)

    def wf_norm(self, wave_fn):
        """
        Returns the norm of a wave function.

        Parameters
        ----------
        wave_fn : array
            Length-N array of the wavefunction in the position representation
        """
        assert wave_fn.shape == self.x.shape
        return np.sqrt((abs(wave_fn) ** 2).sum() * 2 * np.pi / self.dx)

    def solve(self, dt, Nsteps=1, eps=1e-3, max_iter=1000):
        """
        Propagate the Schrodinger equation forward in imaginary
        time to find the ground state.

        Parameters
        ----------
        dt : float
            The small time interval over which to integrate
        Nsteps : float, optional
            The number of intervals to compute (default = 1)
        eps : float
            The criterion for convergence applied to the norm (default = 1e-3)
        max_iter : float
            Maximum number of iterations (default = 1000)
        """
        eps = abs(eps)
        assert eps > 0
        t0 = self.t
        old_psi = self.psi_x
        d_psi = 2 * eps
        num_iter = 0
        while (d_psi > eps) and (num_iter <= max_iter):
            num_iter += 1
            self.time_step(-1j * dt, Nsteps)
            d_psi = self.wf_norm(self.psi_x - old_psi)
            old_psi = 1. * self.psi_x
        self.t = t0

    def time_step(self, dt, Nsteps=1):
        """
        Perform a series of time-steps via the time-dependent Schrodinger
        Equation.

        Parameters
        ----------
        dt : float
            The small time interval over which to integrate
        Nsteps : float, optional
            The number of intervals to compute.  The total change in time at
            the end of this method will be dt * Nsteps (default = 1)
        """
        assert Nsteps >= 0
        self.dt = dt
        if Nsteps > 0:
            self.psi_mod_x *= self.x_evolve_half
            for num_iter in xrange(Nsteps - 1):
                self.compute_p_from_x()
                self.psi_mod_p *= self.p_evolve
                self.compute_x_from_p()
                self.psi_mod_x *= self.x_evolve
            self.compute_p_from_x()
            self.psi_mod_p *= self.p_evolve
            self.compute_x_from_p()
            self.psi_mod_x *= self.x_evolve_half
            self.compute_p_from_x()
            self.psi_mod_x /= self.norm
            self.compute_p_from_x()
            self.t += dt * Nsteps
