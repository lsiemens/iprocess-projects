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
        
        self.fig = None
        self.ax = None
        self.lines = None
        self.anim = None
        self.sparse = 1

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
        
    def setup_animation(self):
        self.fig, self.ax = pyplot.subplots()
        self.lines = []
        for i in xrange(self.n):
            line, = self.ax.plot([],[])
            self.lines.append(line)
#        self.ax.set_xlim([-self.n*self.l*1.1, self.n*self.l*1.1])
#        self.ax.set_ylim([-self.n*self.l*1.1, self.n*self.l*1.1])
        self.ax.set_xlim([-0.0004, 0.0004])
        self.ax.set_ylim([-1.1, 0.1])
        
    def animate(self, i):
        for j in xrange(self.n):
            self.lines[j].set_data([data.x[i*self.sparse][j], data.x[i*self.sparse][j+1]], [data.y[i*self.sparse][j], data.y[i*self.sparse][j+1]])
        return tuple(self.lines)

    def show_animation(self, animate, sparse = 1):
        self.sparse = sparse
        self.anim = animation.FuncAnimation(self.fig, animate, frames=int(data.sets*data.set_size/float(self.sparse)), interval=20)
        pyplot.show()
        
data = npendulum()
data.get_data("output.dat")

data.setup_animation()

def animate(i):
    return data.animate(i)

data.show_animation(animate, sparse = 50)
