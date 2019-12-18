"""
MARkov STAts: check for symbol dependence using Markov chains
"""

"""
STYLE GUIDE

My numpydoc description of a kind
of very exhautive numpydoc format docstring.

Parameters
----------
first : array_like
    the 1st param name `first`
second :
    the 2nd param
third : {'value', 'other'}, optional
    the 3rd param, by default 'value'

Returns
-------
string
    a value in a string

Raises
------
KeyError
    when a key error
OtherError
    when an other error
"""

import sys
import numpy
import scipy
import scipy.stats
import itertools

class MarSta:
    def __init__(self, raw_data, order=5):
        """
        Parameters
        ----------
        self : self
        raw_data : iterable
            sequence of symbols to be analyzed
        order : integer=5
            order of the Markov chain, the analasis considers the previous (order) symbols.
        """
        self.order = order
        self._raw_data = raw_data
        self.symbols = list(set(self._raw_data))
        self.data = {}

        self.train()

    def train(self):
        self.data = {}
        for i in range(len(self._raw_data) - self.order):

            if len(self._raw_data) - self.order > 100:
                if i % ((len(self._raw_data) - self.order)//100) == 0:
                    sys.stdout.write("\r%d%%" % int(100*i/(len(self._raw_data) - self.order)))
                    sys.stdout.flush()

            substring = self._raw_data[i:i + self.order + 1]
            for j in range(self.order + 1):
                key = tuple(substring[j:])
                if key in self.data:
                    self.data[key] += 1
                else:
                    self.data[key] = 1
        print(", Training complete.\n")

    # TODO add allUniform check

    def allIndependent(self):
        data = []
        for i in range(1, self.order + 1):
            for ngram in itertools.product(self.symbols, repeat=i):
                if ngram in self.data:
#                    try:
                        data.append((ngram, self.isIndependent(ngram)))
#                    except ValueError:
#                        print("inssuficient data for ngram ", ngram)
        return data

    def isUniform(self, ngram=()):
        distribution = numpy.array([1.0/len(self.symbols)]*len(self.symbols))
        return self.isDist(distribution, ngram)

    def isIndependent(self, ngram):
        ngram = tuple(ngram)
        if ngram is ():
            raise ValueError("Error: cannot test against itself")
        options = self.getOptions()
        values = [option[1] for option in options]
        if any([value == 0 for value in values]):
            raise ValueError("Error: one or more symbols occured zero times.")
        distribution = numpy.array(values)/sum(values)
        return self.isDist(distribution, ngram)

    def isDist(self, distribution, ngram=()):
        ngram = tuple(ngram)
        options = self.getOptions(ngram)
        values = [option[1] for option in options]
#        if any([value == 0 for value in values]):
#            raise ValueError("Error: one or more symbols occured zero times.")

        values = numpy.array(values)
        mask = numpy.array(values)
        mask[mask<5] = 0
        mask[mask>=5] = 1
        mask = mask.astype(numpy.bool)
        values = values[mask]
        num = sum(values)
        model = numpy.array(distribution)
        model = model[mask]
        model = num*model/sum(model)
        print(values)
        print(model)
#        if any([value < 5 for value in model]):
#            raise ValueError("Error: insufficent samples")
        return scipy.stats.chisquare(f_obs=values, f_exp=model)

    def getOptions(self, ngram=()):
        """
        Parameters
        ----------
        self : self
        ngram : tuple=()
            The previous len(ngram) symbols to consider when finding observed options.

        Returns
        -------
        dict
            A dictonary of all options
        """
        ngram = tuple(ngram)
        options = []
        for symbol in self.symbols:
            key = ngram + (symbol,)
            if key in self.data:
                options.append((symbol,self.data[key]))
            else:
                options.append((symbol,0))
        return options
