####
#
# Copyright (c) 2015, Luke Siemens
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright 
# notice, this list of conditions and the following disclaimer in the 
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its 
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
####

"""
analytic soutions to the schrodingers equation

AUTHOR: Luke Siemens
"""

import cmath
import numpy as np
import numpy.polynomial.hermite as hermite
from matplotlib import pyplot

class analytic_solution:
    def __init__(self, x, m=1.0, dt=0.01, L=None):
        self.x = np.asarray(x)
        self.dx = self.x[1]-self.x[0]

        if L == None:
            self.L = float(self.x[-1] - self.x[0])
        else:
            self.L = float(L)

        self.t = 0.0
        self.dt = dt
        self.m = float(m)
        self.N = len(x)
        self.V_x = np.zeros(self.x.shape)
        self.Cns = np.array([], dtype=complex)
        self._cache = [(None,0)] #empty first element

    def clear_cache(self):
        self._cache = [(None,0)]
    
    def get_psi_n(self, n):
        raise NotImplementedError("get_psi_n not implimented")

    def get_energy_n(self, n):
        raise NotImplementedError("get_energy_n not implimented")

    def time_step(self):
        self.t += self.dt
        
    def get_axis(self):
        return self.x

    def get_psi(self):
        psi_x = np.zeros(self.x.shape, dtype=complex)
        for n, Cn in enumerate(self.Cns):
            if n != 0:
                if len(self._cache) - 1 >= n:
                    psi_x = psi_x + Cn*self._cache[n][0]*np.exp(-1j*self._cache[n][1]*self.t)
                else:
                    psi_n = self.get_psi_n(n)
                    energy_n = self.get_energy_n(n)
                    self._cache.append((psi_n, energy_n))
                    psi_x = psi_x + Cn*psi_n*np.exp(-1j*energy_n*self.t)
        return psi_x
        
    def add_eigenstate(self, n, Cn):
        try:
            if n + 1 > len(self.Cns):
                new_Cns = np.zeros((n + 1,))
                new_Cns[:len(self.Cns)] = self.Cns
                self.Cns = new_Cns
            self.Cns[n] = Cn
        except TypeError:
            if np.max(n) + 1 > len(self.Cns):
                new_Cns = np.zeros((np.max(n) + 1,), dtype=complex)
                new_Cns[:len(self.Cns)] = self.Cns
                self.Cns = new_Cns
            for i in xrange(len(Cn)):
                self.Cns[n[i]] = Cn[i]
        self.Cns *= 1/np.sqrt(np.sum(np.conj(self.Cns)*self.Cns))

    def eigenbasis(self, n_max, psi_x):
        assert psi_x.shape == self.x.shape
        self.Cns = np.zeros((n_max,), dtype=complex)
        for n in xrange(n_max):
            if n!=0:
                integral = self.dx*np.sum(np.multiply(np.conj(self.get_psi_n(n)), psi_x))
                self.Cns[n] = integral
        self.Cns *= 1/np.sqrt(np.sum(np.conj(self.Cns)*self.Cns))

class inf_square_well(analytic_solution):
    def get_psi_n(self, n):
        if n%2 == 1:
            psi = np.sqrt(2/(self.L))*np.cos(n*np.pi*self.x/self.L)
        else:
            psi = np.sqrt(2/(self.L))*np.sin(n*np.pi*self.x/self.L)
        psi[self.x > self.L/2.0] = 0
        psi[self.x < -self.L/2.0] = 0
        return psi

    def get_energy_n(self, n):
        return (n*np.pi/self.L)**2/(2*self.m)

class harmonic_well(analytic_solution):
    def __init__(self, x, k=1.0, m=1.0, dt=0.01, L=None):
        analytic_solution.__init__(self, x, m, dt, L)
        self.k = k
        self.omega = np.sqrt(self.k/self.m)

    def get_psi_n(self, n):
        psi = np.sqrt(1/float(np.math.factorial(n-1)*2**(n-1)))*(self.m*self.omega/np.pi)**(0.25)*np.exp(-self.m*self.omega*self.x**2/2)*hermite.hermval(np.sqrt(self.m*self.omega)*self.x,[0]*(n-1)+[1.0])
        psi[self.x > self.L/2.0] = 0
        psi[self.x < -self.L/2.0] = 0
        return psi

    def get_energy_n(self, n):
        return self.omega*(n-0.5)
