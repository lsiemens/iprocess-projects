import sys
import numpy

class MarkovTest:
    def __init__(self, raw_data, order=5):
        self.order = order
        self.raw_data = raw_data
        self._data = {}
        self.symbols = list(set(raw_data))

        self.calulate_chain()

    def calculate_chain(self):
        self._data = {}
        for i in range(len(self.raw_data) - self.order):

            if i % ((len(self.raw_data) - self.order)//100) == 0:
                sys.stdout.write("\r%d%%" % int(100*i/(len(self.raw_data) - self.order)))
                sys.stdout.flush()

            subset = self.raw_data[i:i + self.order + 1]
            for j in range(self.order + 1):
                key = tuple(subset[j:])
                if key in self._data:
                    self._data[key] += 1
                else:
                    self._data[key] = 1
        print(", complete.\n")

    def getOptions(self, ngram=()):
        options = []
        for symbol in self.symbols:
            key = ngram + tuple([symbols])
            if key in self._data:
                options.append((symbol, self._data[key]))
            else:
                options.append((symbol, 0))
        return options
