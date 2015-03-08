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

import numpy
from matplotlib import animation
from matplotlib import pyplot

class plotting:
    def __init__(self, n_x=1, n_y=1):
        self.n_x = n_x
        self.n_y = n_y

        self.lines = []
        self.animator = None
        
        self.fig = pyplot.figure()
        self.axis = [[None for j in xrange(n_x)] for i in xrange(n_y)]
        for i in xrange(n_y):
            for j in xrange(n_x):
                self.axis[i][j] = self.fig.add_subplot(n_y, n_x, (j + 1)+n_x*i)
                self.axis[i][j].set_title(str((j+1)+n_x*i))
    
    def add_line(self, n_x, n_y, get_x, get_y, time_step, form="prob"):
        axis_line, = self.axis[n_y - 1][n_x - 1].plot(get_x(), get_y())
        self.lines.append([n_x, n_y, get_x, get_y, time_step, axis_line, form])

    def _animate_init(self):
        data = ()
        for line in self.lines:
            line[5].set_data([], [])
            if line[4] != None:
                line[4]()
            data = data + (line[5],)
        return data

    def _animate_plot(self, i):
        data = ()
        for line in self.lines:
            if line[6] == "prob":
                line[5].set_data(line[2](), numpy.conj(line[3]())*line[3]())
            elif line[6] == "real":
                line[5].set_data(line[2](), numpy.real(line[3]()))
            elif line[6] == "imag":
                line[5].set_data(line[2](), numpy.imag(line[3]()))
            if line[4] != None:
                line[4]()
            data = data + (line[5],)
        return data

    def animate(self):
        self.animator = animation.FuncAnimation(self.fig, self._animate_plot, init_func=self._animate_init, frames=100, interval=30, blit=True)
        pyplot.show()

    def plot(self):
        for line in self.lines:
            n_x = line[0]
            n_y = line[1]
            self.axis[n_y - 1][n_x - 1].plot(line[2](), line[3]())
        pyplot.show()

import analytic
import schrodinger

dt = 10000
N = 2**16
M = 2**11
dx = 1.0
x = dx*(numpy.arange(N) - 0.5*N)

x_lim = dx*M

V_x = numpy.zeros(x.shape)
V_x[x < -x_lim] = 1E-5
V_x[x > x_lim] = 1E-5
V_x[x < -x_lim] = 1E100
V_x[x > x_lim] = 1E100

analytic = analytic.inf_square_well(x=x, m=1, dt=dt, L=2*x_lim)
analytic.add_eigenstate([1, 2], [1, 1])

psi_x0 = analytic.get_psi()

numeric = schrodinger.Schrodinger(x=x, psi_x0=psi_x0, V_x=V_x, m=1)

def time_step():
    numeric.time_step(dt)

def potential():
    return V_x
    
def pax():
    return numeric.p

plot = plotting(3,2)
plot.add_line(1, 1, analytic.get_axis, analytic.get_psi, analytic.time_step, 'prob')
plot.add_line(1, 2, analytic.get_axis, analytic.get_psi, None, 'real')
plot.add_line(1, 2, analytic.get_axis, analytic.get_psi, None, 'imag')

plot.add_line(2, 1, analytic.get_axis, numeric._get_psi_x, time_step, 'prob')
plot.add_line(2, 2, analytic.get_axis, numeric._get_psi_x, None, 'real')
plot.add_line(2, 2, analytic.get_axis, numeric._get_psi_x, None, 'imag')

plot.add_line(3, 1, pax, numeric._get_psi_p, time_step, 'prob')
plot.add_line(3, 2, pax, numeric._get_psi_p, None, 'real')
plot.add_line(3, 2, pax, numeric._get_psi_p, None, 'imag')

#plot.add_line(2, 1, analytic.get_axis, potential, None, 'real')
plot.animate()
