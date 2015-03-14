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
        self.find_energy = False
        self.sets = 0
        self.set_size = 0
        self.n = 0
        self.l = 0
        self.g = 0
        self.dt = 0
        self.theta = None
        self.theta_d = None
        self.theta_dd = None
        self.energy = None
        self.x = None
        self.y = None
        
        self.fig = None
        self.ax = None
        self.lines = None
        self.time_sig = None
        self.energy_sig = None
        self.anim = None
        self.sparse = 1
        
        self.time_format = "{:1.4E}s"
        self.energy_format = "{:1.4E}j"

    def get_data(self, fname):
        file_in = pyio.pyio(fname)
        file_in.load()
        self.format = file_in.format
        self.find_energy = file_in.find_energy
        self.sets = file_in.sets
        self.set_size = file_in.set_size
        self.n = file_in.n
        self.l = file_in.l
        self.g = file_in.g
        self.dt = file_in.dt
        data = file_in.data
        
        if self.find_energy:
            self.energy = numpy.array(file_in.energy)
            
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
        
    def setup_animation(self, window_min = 1):
        self.fig, self.ax = pyplot.subplots()
        self.lines = []
        x_max = numpy.max(self.x)
        x_min = numpy.min(self.x)
        y_max = numpy.max(self.y)
        y_min = numpy.min(self.y)

        x_avg = (x_max + x_min)/2.0
        y_avg = (y_max + y_min)/2.0

        x_width = (x_max - x_min)/2.0
        y_width = (y_max - y_min)/2.0
        
        if x_width < window_min/2.0:
            x_width = window_min/2.0

        if y_width < window_min/2.0:
            y_width = window_min/2.0
        
        if numpy.isnan(x_avg):
            x_avg = 0.0
            x_width = self.n*self.l

        if numpy.isnan(y_avg):
            y_avg = 0.0
            y_width = self.n*self.l
            
        x_width = x_width*1.1
        y_width = y_width*1.1
        
        for i in xrange(self.n):
            line, = self.ax.plot([],[])
            self.lines.append(line)

        self.time_sig = self.ax.text(0.05, 0.95, "time: ", transform=self.ax.transAxes)
        
        if self.find_energy:
            self.energy_sig = self.ax.text(0.05, 0.85, "Energy: ", transform=self.ax.transAxes)

        self.ax.set_xlim([x_avg-x_width, x_avg+x_width])
        self.ax.set_ylim([y_avg-y_width, y_avg+y_width])
        
    def animate(self, i):
        for j in xrange(self.n):
            self.lines[j].set_data([data.x[i*self.sparse][j], data.x[i*self.sparse][j+1]], [data.y[i*self.sparse][j], data.y[i*self.sparse][j+1]])
        self.time_sig.set_text("time: " + self.time_format.format(self.dt*i))
        if self.find_energy:
            self.energy_sig.set_text("Energy: " + self.energy_format.format(self.energy[i]))
            return tuple(self.lines)
        else:
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
