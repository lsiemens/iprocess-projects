#!/usr/bin/env ipython3

from matplotlib import pyplot
import numpy
import os
import text_gen

data = text_gen.markov_chain()
data.load("sav.dat")

dist_max = None
distribution = {}
for key, value in data._transitions.items():
    if dist_max is not None:
        if dist_max < value:
            dist_max = value
    else:
        dist_max = value
        
    if value in distribution:
        distribution[value] += 1
    else:
        distribution[value] = 1

data = numpy.zeros(shape=(dist_max + 1))

print(distribution)

for key, value in distribution.items():
    data[key] = value

pyplot.plot(data)
pyplot.show()
