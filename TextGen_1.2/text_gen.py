"""

text generator 1.02

Text is generated using a markov chain.

using sparse data storage.

"""

from itertools import izip
import numpy
import random
from time import time
import cPickle
import sdata

class markov_chain:
    def __init__(self):
#        self.valid_chars = "cepwoftdn, =1234567890\n"
        self.valid_chars = "abcdefghijklmnopqrstuvwxyz. "
        self.order = 0
        self.transitions = sdata.sdata()
    
    
    def calculate_transitions(self, file_name, order, calculate_transitions=True):
        self.order = order
        self.file_name = file_name
        if isinstance(self.file_name, basestring):
            self.file_name = [self.file_name]

        text = ""
        for fname in self.file_name:
            with open(fname, 'r') as txt_file:
#                tmp = "".join(char for char in txt_file.read().lower() if char in self.valid_chars + "\n")
#                text += tmp.replace("\n", " ")
                text += "".join(char for char in txt_file.read().lower() if char in self.valid_chars)
        self.text = text
        self.transitions = sdata.sdata((len(self.valid_chars), )*(order + 1))
        
        if calculate_transitions:
            substring = []
            print "\ngenerating the transition space.\n"
            for i in xrange(len(text) - self.order):
                if (i) % ((len(text) - self.order)/10) == 0:
                    percentage = float(i*100)/(len(text) - self.order)
                    print ("%.1f" % percentage) + " %"
                substring = text[i:i + self.order + 1]
                indices = tuple([self.valid_chars.index(char) for char in substring])
                self.transitions[indices] += 1
            print "\nnormalizing the transition space.\n"
            self.transitions._normalize(self.valid_chars)

    def generate_string(self, length, seed=""):
        if seed == "":
            seed = self.text[:self.order]
        if len(seed) < self.order:
            raise ValueError("the length of the seed must be equal to the order")
        text = seed
        while (len(text) < length):
            if self.order != 0:
                context = text[-self.order:]
                indices = tuple([self.valid_chars.index(char) for char in context])
                text += self.choose_random(self.valid_chars, self.transitions[indices])
            else:
                text += self.choose_random(self.valid_chars, self.transitions)
        with open("out.txt", 'w') as txt_file:
            txt_file.write(text)
        print text
        
    def choose_random(self, items, probabilities):
        x = random.uniform(0, 1)
        total_probability = 0.0
        for item, probability in zip(items, probabilities):
            total_probability += probability
            if x < total_probability:
                return item
        return items[-1]
    
    def save(self, file_name):
        with open(file_name, 'wb') as data_file:
            print "saving state."
            pickler = cPickle.Pickler(data_file, -1)
            pickler.dump(self.order)
            pickler.dump(self.valid_chars)
#            pickler.dump(self.transitions)
        self.transitions.save(file_name)

    def load(self, file_name):
        with open(file_name, 'rb') as data_file:
            print "loading state."
            unpickler = cPickle.Unpickler(data_file)
            self.order = unpickler.load()
            self.valid_chars = unpickler.load()
#            self.transitions = unpickler.load()
            self.transitions.load(unpickler)
    
    def _normalize(self, array, first=True):
        try:
            len(array[0])
            if first:
                print "0.0 %"
            for i, sub_array in enumerate(array):
               self._normalize(sub_array, first=False)
               if first:
                   percentage = 100*(i + 1)/float(len(self.valid_chars))
                   print ("%.1f"%percentage) + " %"
            if first:
                print "\n" 
        except TypeError:
            sum = numpy.sum(array)
            if sum != 0:
                array[:] /= sum

dir = "raw_text//HGWells//"
#dir = "raw_text//"
ext = ".txt"
order = 8
transitions_file = "transition_hgwells_order_"+str(order)
calculate = True

names = open(dir+"f_names.txt", 'r').read().split('\n')
#names = ["portal"]
text_gen = markov_chain()

if calculate:
    text_gen.calculate_transitions([dir+name+ext for name in names], order)
#    text_gen.save(dir + transitions_file + ".data")
#else:
#    text_gen.calculate_transitions([dir+name+ext for name in names], 0, False)
#    text_gen.load(dir + transitions_file + ".data")
    text_gen.generate_string(2000)
