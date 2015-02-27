"""
General Numerical Solver for the 1D Time-Dependent Schrodinger Equation.

Authors:
- Jake Vanderplas <vanderplas@astro.washington.edu>
- Andre Xuereb (imaginary time propagation, normalized wavefunction
- Luke Siemens <luke@lsiemens.com> (improved imaginary time propogation)

For a theoretical description of the algorithms, please see
http://jakevdp.github.com/blog/2012/09/05/quantum-python/
and http://lsiemens.com/psipy-solving-schrodingers-equation/

For the original version of this code please see
http://github.com/jakevdp/pySchrodinger/

This version of schrodinger.py has been modifi to produce results in
momentum-space not wavenumber-space, and the method for finding the
ground state has been rewriten to numericaly calculate successive states
of the hamiltonian operator.

License: BSD style
Please feel free to use and modify this, but keep the above information.
"""

import numpy as np
import numpy.ma as ma
from scipy import fftpack

import matplotlib.pyplot as pyplot

class Schrodinger(object):
    """
    Class which implements a numerical solution of the time-dependent
    Schrodinger equation for an arbitrary potential
    """
    def __init__(self, x, psi_x0, V_x, m=1):
        """
        Parameters
        ----------
        x : array_like, float
            Length-N array of evenly spaced spatial coordinates
        psi_x0 : array_like, complex
            Length-N array of the initial wave function at time t0
        V_x : array_like, float
            Length-N array giving the potential at each x
        m : float
            Particle mass (default = 1)
        """
        # Validation of array inputs
        self.x, psi_x0, self.V_x = map(np.asarray, (x, psi_x0, V_x))
        N = self.x.size
        assert self.x.shape == (N,)
        assert psi_x0.shape == (N,)
        assert self.V_x.shape == (N,)

        # Validate and set internal parameters
        assert m > 0
        self.m = m
        self.t = 0.0
        self.dt_ = None
        self.N = len(x)
        self.dx = self.x[1] - self.x[0]
        self.dp = 2 * np.pi / (self.N * self.dx)
        self._near_zero = 1e-10 #values below this are considered near zero

        # Set momentum scale
        self.p0 = -np.pi/self.dx
        self.p = self.p0 + self.dp * np.arange(self.N)
        
        self.psi_x = psi_x0
        self.compute_p_from_x()

        # Variables which hold steps in evolution
        self.x_evolve_half = None
        self.x_evolve = None
        self.p_evolve = None

    def _set_psi_x(self, psi_x, normalize=True):
        assert psi_x.shape == self.x.shape
        self.psi_mod_x = (psi_x * np.exp(-1j * self.p[0] * self.x)
                          * self.dx / np.sqrt(2 * np.pi))
        if normalize:
            self.normalize()
        self.compute_p_from_x()

    def _get_psi_x(self):
        return (self.psi_mod_x * np.exp(1j * self.p[0] * self.x)
                * np.sqrt(2 * np.pi) / self.dx)

    def _set_psi_p(self, psi_p, normalize=True):
        assert psi_p.shape == self.x.shape
        self.psi_mod_p = psi_p * np.exp(1j * self.x[0] * self.dp
                                        * np.arange(self.N))
        self.compute_x_from_p()
        if normalize:
            self.normalize()

    def _get_psi_p(self):
        return self.psi_mod_p * np.exp(-1j * self.x[0] * self.dp
                                        * np.arange(self.N))

    def _get_dt(self):
        return self.dt_

    def _set_dt(self, dt):
        assert dt != 0
        if dt != self.dt_:
            self.dt_ = dt
            self.x_evolve_half = np.exp(-0.5 * 1j * self.V_x * self.dt)
            self.x_evolve = self.x_evolve_half ** 2
            self.p_evolve = np.exp(-0.5 * 1j * (self.p ** 2) * self.dt
                                    / (self.m))

    def normalize(self):
        self.psi_mod_x *= self.wf_norm(self.psi_x)
        self.compute_p_from_x()

    psi_x = property(_get_psi_x, _set_psi_x)
    psi_p = property(_get_psi_p, _set_psi_p)
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
        return 1/np.sqrt(self.dx*np.sum(np.multiply(np.conj(wave_fn), wave_fn)))

    def hamiltonian_eigenstate(self, dt, eigenstates=[], Nsteps=1, eps=1e-3, max_iter=1000):
        """
        Propagate the Schrodinger equation in imaginary
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
        eigenstates = np.array(eigenstates)
        t0 = self.t
        psi_x0 = np.copy(self.psi_x)
        
        old_psi = ma.masked_less(np.multiply(np.conj(self.psi_x), self.psi_x), self._near_zero)
        decay_variance = 2 * eps #the 
        num_iter = 0
        while (decay_variance > eps):
            if num_iter > max_iter:
                self.t = t0
                self.psi_x = psi_x0
                self.compute_p_from_x()
                raise RuntimeError("faild to converge to an eigenstate after " + str(num_iter - 1) + " iterations.")
            num_iter += 1
            
            self.time_step(-1j * dt, Nsteps, normalize = False)
            if len(eigenstates) > 0:
                for i, eigenstate in enumerate(eigenstates):
                    Cn = np.sum(np.multiply(self.psi_x, eigenstate))*self.dx
                    self._set_psi_x(self.psi_x - Cn*eigenstate, normalize=False)
                self.compute_p_from_x()
            
            mask_psi_x = ma.masked_less(np.multiply(np.conj(self.psi_x), self.psi_x), self._near_zero)
            decay = np.real(mask_psi_x/old_psi) #both values should be real bu just in case force it to be real
            decay_variance = np.var(decay)
            print decay_variance
            #manualy normalize
            self.normalize()
            
            old_psi = mask_psi_x
            pyplot.title("decay")
            pyplot.plot(self.x, decay, c='k')
            pyplot.show()
        #ma.log is used rather than np.log so it can handle the zeros in decay
        energy = -(1.0/(2*dt*Nsteps))*np.mean(ma.log(decay))
        denergy = (1.0/(2*dt*Nsteps))*np.std(ma.log(decay))
        print np.std((-1.0/(2*self.m))*(ma.masked_less(np.diff(self.psi_x,2)/(self.dx**2), self._near_zero)/ma.masked_less(self.psi_x[1:-1], self._near_zero)) + self.V_x[1:-1])

#        energy = 0
#        denergy = 0
        eigenstate = np.copy(self.psi_x) 
        self.t = t0
        self.psi_x = psi_x0
        self.compute_p_from_x()
        return eigenstate, (energy, denergy)
        

    def time_step(self, dt, Nsteps=1, normalize = True):
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
            if normalize:
                self.normalize()
            self.compute_p_from_x()
            self.t += dt * Nsteps
