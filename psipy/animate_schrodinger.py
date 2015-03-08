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


def square_barrier(x, width, height):
    return height * (theta(x) - theta(x - width))

######################################################################
# Create the animation

unit_sys = units.units(10**-12, mode = "abs")
_T = unit_sys.get_T()
_T.set_format("{:1.3e}")

# specify time steps and duration
dt = 0.00001
N_steps = 50
t_max = 1000000
ylim = (0, 1.0)
frames = 100

# specify constants
m = 1.0      # particle mass

# specify range in x coordinate
fac = 1
N = 2 ** (11 + fac)
dx = 0.005*2**(-fac)
x = dx * (np.arange(N) - 0.5 * N)
print x

# specify potential
V_x = square_barrier(x, 1.5, 0)#.5)
V_x[x < -5] = 1E3
V_x[x > 5] = 1E3

psi_x0 = gauss_x(x, 0.1, 0, 10.0)
#psi_x0 = np.cos(np.pi*x/(2*98))
#psi_x0 =psi_x0 + np.cos(3*np.pi*x/(2.0*98))

# define the Schrodinger object which performs the calculations
S = Schrodinger(x=x,
                psi_x0=psi_x0,
                V_x=V_x,
                m=m)

######################################################################
# Set up plot
fig = plt.figure()

# plotting limits
xlim = (-8, 8)
print S.p
plim = (-100, 100)

# top axes show the x-space data
ax1 = fig.add_subplot(211, xlim=xlim, ylim=ylim)
psi_x_line, = ax1.plot(S.x, S.psi_x, c='r', label=r'$|\psi(x)|$')
V_x_line, = ax1.plot(S.x, S.psi_x, c='k', label=r'$V(x)$')
center_line = ax1.axvline(0, c='k', ls=':', label=r"$x_0 + v_0t$")

time = ax1.text(0,0,"")
ax1.legend(prop=dict(size=12))
ax1.set_xlabel('$x$')
ax1.set_ylabel(r'$|\psi(x)|$')

# bottom axes show the k-space data
ax2 = fig.add_subplot(212, xlim=plim, ylim=(-1, 5.5))
psi_p_line, = ax2.plot([], [], c='r', label=r'$|\psi(p)|$')

ax2.legend(prop=dict(size=12))
ax2.set_xlabel('$p$')
ax2.set_ylabel(r'$|\psi(p)|$')

######################################################################
# Functions to Animate the plot
def init():
    psi_x_line.set_data([], [])
    V_x_line.set_data([], [])

    psi_p_line.set_data([], [])
    time.set_text("")
    return (psi_x_line, V_x_line, psi_p_line, time)

def animate(i):
    S.time_step(dt, N_steps)
    psi_x_line.set_data(S.x, np.real(np.conj(S.psi_x)*S.psi_x))
    V_x_line.set_data(S.x, S.V_x)

    psi_p_line.set_data(S.p, abs(S.psi_p))
    time.set_text("t = " + str(abs(S.t)*_T))
    return (psi_x_line, V_x_line, psi_p_line, time)

# call the animator.
# blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=frames, interval=30, blit=True)


# uncomment the following line to save the video in mp4 format.  This
# requires either mencoder or ffmpeg to be installed on your system
#anim.save('schrodinger_barrier.mp4', fps=15,
#          extra_args=['-vcodec', 'libx264'])

plt.show()
