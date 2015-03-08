####
#
# Copyright (c) 2015, Jake Vanderplas, Luke Siemens
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
Solve and animate the Schrodinger equation

First presented at http://jakevdp.github.com/blog/2012/09/05/quantum-python/

Authors:
- Jake Vanderplas <vanderplas@astro.washington.edu>
- Luke Siemens (small modifications, switching from k to p-space)

License: BSD
"""

import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
from schrodinger import Schrodinger
import units

######################################################################
# Helper functions for gaussian wave-packets
def gauss_x(x, a, x0, k0):
    """
    a gaussian wave packet of width a, centered at x0, with momentum k0
    """
    return ((a * np.sqrt(np.pi)) ** (-0.5)
            * np.exp(-0.5 * ((x - x0) * 1. / a) ** 2 + 1j * x * k0))


def gauss_p(p, a, x0, p0):
    """
    analytical fourier transform of gauss_x(x), above
    """
    return ((a / np.sqrt(np.pi)) ** 0.5
            * np.exp(-0.5 * (a * (p - p0)) ** 2 - 1j * (p - p0) * x0))

######################################################################
# Utility functions for running the animation
def theta(x):
    """
    theta function :
      returns 0 if x<=0, and 1 if x>0
    """
    x = np.asarray(x)
    y = np.zeros(x.shape)
    y[x > 0] = 1.0
    return y


def square_barrier(x, width, height, x0=0.0):
    return height * (theta(x-x0) - theta(x - width - x0))

######################################################################
# Create the animation
unit_sys = units.units(10**-12, mode="abs")
_T = unit_sys.get_T()
_T.set_format("{:1.3e}")
_E = unit_sys.get_E()

# specify time steps and duration
dt = 1.0
N_steps = 50
t_max = 120
frames = int(t_max / float(N_steps * dt))

# specify constants
m = 1      # particle mass

# specify range in x coordinate
N = 2 ** 11
dx = 0.1
x = dx * (np.arange(N) - 0.5 * N)

# specify potential
V0 = 1.5
L = 1.0 / np.sqrt(2 * m * V0)
a = 3 * L
x0 = -60 * L
V_x = square_barrier(x, a, 0)
V_x[x < -98] = 1E30
V_x[x > 98] = 1E30
V_x[x < -70] = 1E-6
V_x[x > 70] = 1E-6

# specify initial momentum and quantities derived from it
p0 = np.sqrt(2 * m * 0.2 * V0)
dp2 = p0 * p0 * 1. / 80
d = 1.0 / np.sqrt(2 * dp2)

v0 = p0 / m
#psi_x0 = gauss_x(x, d, x0, p0)
psi_x0 = np.cos(5*np.pi*x/(2*98))
psi_x0 = psi_x0 + np.cos(21*np.pi*x/(2*98))

# define the Schrodinger object which performs the calculations
S = Schrodinger(x=x,
                psi_x0=psi_x0,
                V_x=V_x,
                m=5*m)

n = 2
state = []
energy = []
try:
    for i in  xrange(n):
        s, E = S.hamiltonian_eigenstate(dt, state, Nsteps=100, eps=1e-9, max_iter = 10000)
        plt.plot(S.x, np.real(s))
        state.append(s)
        energy.append(E)
    for E in energy:
        print E[0]*_E, E[1]*_E 
    plt.show()
except RuntimeError:
    print "not convergent"

######################################################################
# Set up plot
fig = plt.figure()

# plotting limits
xlim = (-100, 100)
plim = (-5, 5)

# top axes show the x-space data
ymin = 0
ymax = V0
ax1 = fig.add_subplot(111, xlim=xlim,
                      ylim=(ymin - 0.2 * (ymax - ymin),
                            ymax + 0.2 * (ymax - ymin)))
psi_x_line, = ax1.plot(S.x, 4*abs(S.psi_x), c='r', label=r'$|\psi(x)|$')
V_x_line, = ax1.plot(S.x, S.V_x, c='k', label=r'$V(x)$')

ax1.set_title("t = " + str(abs(S.t)*_T))
ax1.legend(prop=dict(size=12))
ax1.set_xlabel('$x$')
ax1.set_ylabel(r'$|\psi(x)|$')
plt.show()
