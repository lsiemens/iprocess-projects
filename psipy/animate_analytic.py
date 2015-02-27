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
import analytic
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
_E = unit_sys.get_E()

# specify time steps and duration
dt = 0.00041
N_steps = 100
t_max = 120
ylim = (-2.0, 2.0)
frames = int(100)

# specify constants
hbar = 1.0   # planck's constant
m = 2.0      # particle mass

# specify range in x coordinate
N = 2 ** 13
dx = 6.0 / float(N)
x = dx * (np.arange(N) - 0.5 * N)

# specify potential
V0 = 1.5
L = x[-1]-x[0]
a = 3 * L
x0 = -60 * L
V_x = square_barrier(x, a, 0.5)
#V_x[x < -0.49] = 1E6
#V_x[x > 0.49] = 1E6

# specify initial momentum and quantities derived from it
p0 = np.sqrt(2 * m * 0.2 * V0)
dp2 = p0 * p0 * 1. / 80
d = hbar / np.sqrt(2 * dp2)

v0 = p0 / m
psi_x0 = gauss_x(x, 10, -50, 0.5)

# define the Schrodinger object which performs the calculations
S = analytic.harmonic_well(x=x, k=4000.0, m=m, dt=dt, L=L)

#S.eigenbasis(130, np.sqrt(2/float(L))*np.cos(11*np.pi*x/L+0.2))
S.eigenbasis(50, gauss_x(x, 0.15, 0, 0))
#S.eigenbasis(151, square_barrier(x+L/4.0, L/2, 1))
#S.add_eigenstate([1,3,6,7,4], [1.0,2.0,3.5,4.1,5.7])
#S.add_eigenstate([3,4], [1.0,1.0])
#S.add_eigenstate([1], [1.0j])

print 1*_E
print S.get_energy_n(3)*_E
print S.get_energy_n(4)*_E
######################################################################

# Set up plot
fig = plt.figure()

# plotting limits
xlim = (-3, 3)
plim = (-28, 28)

# top axes show the x-space data
ax1 = fig.add_subplot(111, xlim=xlim, ylim=ylim)
psi_x_line, = ax1.plot([], [], c='r', label=r'$|\psi(x)|$')

time = ax1.text(0, 0, "")
ax1.legend(prop=dict(size=12))
ax1.set_xlabel('$x$')
ax1.set_ylabel(r'$|\psi(x)|$')

######################################################################
# Functions to Animate the plot
def init():
    psi_x_line.set_data([], [])
    time.set_text("")
    return (psi_x_line, time)

def animate(i):
    psi = S.animate()
    psi = np.real(np.conj(psi)*psi)
    psi_x_line.set_data(S.x, psi)
    time.set_text("t = " + str(abs(S.t)*_T))
    return (psi_x_line, time)

# call the animator.
# blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=frames, interval=30, blit=True)


# uncomment the following line to save the video in mp4 format.  This
# requires either mencoder or ffmpeg to be installed on your system
#anim.save('schrodinger_barrier.mp4', fps=15,
#          extra_args=['-vcodec', 'libx264'])

plt.show()
