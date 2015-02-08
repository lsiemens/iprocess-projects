"""
A numerical solver for the 2D time dependent schrodinger equation.
"""

import numpy
from scipy import fftpack as fft

class psipy:
    """
    
    Numericaly solve the time dependent shrodinger equation in
    one dimension.

    """
    
    def __init__(self, x, psi_t0, V, k0):
        self._x = x
        self._V = V
        self._k0 = k0
        
        self._N = len(self._x)
        self._t = 0.0
        self._dt = None
        self._dx = self._x[1] - self._x[0]
        self._dk = 2*numpy.pi / (self._N * self._dx)
        self._k = self._k0 + self._dk*numpy.arange(self._N)
        
        self._psi_x = None
        self._psi_k = None
        
        self._x_half_step = None
        self._x_step = None
        self._k_step = None

        self.set_psi_x(psi_t0)
    
    def get_psi_x(self):
        return self._psi_x*numpy.exp(1j*self._k0*self._x)*numpy.sqrt(2*numpy.pi)/self._dx
        
    def set_psi_x(self, psi_x):
        self._psi_x = psi_x*numpy.exp(-1j*self._k0*self._x)*self._dx/numpy.sqrt(2*numpy.pi)
        self._update_k()

    def get_psi_k(self):
        return self._psi_k*numpy.exp(-1j*self._x[0]*self._dk*numpy.arange(self._N))
        
    def set_psi_k(self, psi_k):
        self._psi_k = psi_k*numpy.exp(1j*self._x[0]*self._dk*numpy.arange(self._N))
        self._update_x()
    
    def set_dt(self, dt):
        self._dt = dt
        self._x_half_step = numpy.exp(-0.5j*self._V*self._dt)
        self._x_step = self._x_half_step**2
        self._k_step = numpy.exp(-0.5j*self._dt*self._k**2)
    
    def _update_x(self):
        self._psi_x = fft.ifft(self._psi_k)
    
    def _update_k(self):
        self._psi_k = fft.fft(self._psi_x)
        
    def step(self, steps=1):
        self._psi_x *= self._x_half_step
        for i in xrange(steps - 1):
            self._update_k()
            self._psi_k *= self._k_step
            self._update_x()
            self._psi_x *= self._x_step
        self._update_k()
        self._psi_k *= self._k_step
        self._update_x()
        self._psi_x *= self._x_half_step
        self._update_k()
        self._t += self._dt*steps

    
   
