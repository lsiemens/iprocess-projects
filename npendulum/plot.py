""" 

display the output from npendulum.f90
"""

import numpy
import matplotlib.pyplot as pyplot
import matplotlib.animation as animation
import pyio

class npendulum:
    def __init__(self):
        self.format = 0
        self.sets = 0
        self.set_size = 0
        self.n = 0
        self.l = 0
        self.g = 0
        self.dt = 0
        self.theta = None
        self.theta_d = None
        self.theta_dd = None
        self.x = None
        self.y = None

    def get_data(self, fname):
        file_in = pyio.pyio(fname)
        file_in.load()
        self.format = file_in.format
        self.sets = file_in.sets
        self.set_size = file_in.set_size
        self.n = file_in.n
        self.l = file_in.l
        self.g = file_in.g
        self.dt = file_in.dt
        data = file_in.data
        
        self.theta = []
        self.theta_d = []
        self.theta_dd = []
        for i in xrange(self.sets*self.set_size):
            if self.format == 1:
                self.theta.append(data[i])
            elif self.format == 2:
                self.theta.append(data[self.format*i])
                self.theta_d.append(data[self.format*i+1])
            elif self.format == 3:
                self.theta.append(data[self.format*i])
                self.theta_d.append(data[self.format*i+1])
                self.theta_dd.append(data[self.format*i+2])
            else:
                raise IOError("invalid format.")
        self.theta = numpy.array(self.theta)
        self.theta_d = numpy.array(self.theta_d)
        self.theta_dd = numpy.array(self.theta_dd)
        
        self.x = []
        self.y = []
        for i in xrange(self.sets*self.set_size):
            x = [0]
            y = [0]
            for j in xrange(self.n):
                x.append(x[j]+self.l*numpy.sin(self.theta[i][j]))
                y.append(y[j]-self.l*numpy.cos(self.theta[i][j]))
            self.x.append(x)
            self.y.append(y)
        self.x = numpy.array(self.x)
        self.y = numpy.array(self.y)

data = npendulum()
data.get_data("output.dat")

fig, ax = pyplot.subplots()
line_1, = ax.plot([], [])
line_2, = ax.plot([], [])
ax.set_xlim([-3, 3])
ax.set_ylim([-3, 3])


def animate(i):
    line_1.set_data([0, data.x[i][1]],[0, data.y[i][1]])
    line_2.set_data([data.x[i][1], data.x[i][2]],[data.y[i][1], data.y[i][2]])
    return line_1, line_2

anim = animation.FuncAnimation(fig, animate, frames=data.sets*data.set_size, interval=10)
pyplot.show()
