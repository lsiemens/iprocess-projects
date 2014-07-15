"""

text generator 1.1

Text is generated using a markov chain.

"""

from itertools import izip
import numpy
import random
from time import time

class markov_chain:
    valid_chars = "abcdefghijklmnopqrstuvwxyz. "
    
    def calculate_transitions(self, file_name, order):
        self.order = order
        self.file_name = file_name
        if isinstance(self.file_name, basestring):
            self.file_name = [self.file_name]

        text = ""
        for fname in self.file_name:
            with open(fname, 'r') as txt_file:
                tmp = "".join(char for char in txt_file.read().lower() if char in self.valid_chars + "\n")
                text += tmp.replace("\n", " ")
#                text += "".join(char for char in txt_file.read().lower() if char in self.valid_chars)
        
        self.transitions = numpy.zeros((len(self.valid_chars), )*(order + 1))
            
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
        self._normalize(self.transitions)

    def generate_string(self, length, seed=""):
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
            numpy.save(data_file, [self.transitions, self.order])
    
    def load(self, file_name):
        with open(file_name, 'rb') as data_file:
            tmp = numpy.load(data_file)
            self.transitions = tmp[0]
            self.order = tmp[1]
    
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

dir = "prepared_text//"
ext = ".txt"
names = ["pg1268"]
text_gen = markov_chain()

text_gen.calculate_transitions([dir+name+ext for name in names], 5)
#text_gen.save("transition.data")
#text_gen.load("transition.data")
text_gen.generate_string(1000, "he said")
