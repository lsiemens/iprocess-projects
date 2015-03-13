import pyio
import numpy

format = 1
sets = 1000
set_size = 100
integrator = 1
n = 60
g = 9.82
length = 1
time = 10


l = length / float(n)
dt = time / float(sets*set_size)
theta = [0 for i in xrange(n)]
theta_d = [0 for i in xrange(n)]
theta_dd = [0 for i in xrange(n)]

theta_d[0] = -0.10
theta_d[1] = 0.00
theta_d[2] = 0.10

file_out = pyio.pyio("initalization.dat")
file_out.initalize("initalization.dat", format, sets, set_size, integrator, l, g, dt, theta, theta_d, theta_dd)
