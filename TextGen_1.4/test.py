import markov
import utils

import itertools
import numpy
import random
from matplotlib import pyplot

isfast = False
isfast = True

path = "./Books/*.txt"
texts = utils.open_files(path, markov.Markov.symbols)
files = list(texts.keys())
random.shuffle(files)

text_markov = "\n".join([texts[file] for file in files[:-2]])
text_test = "\n".join([texts[file] for file in files[-2:]])

if isfast:
    print("Fast training-----------------------")
    text_markov = text_markov[10000:110000]
    text_test = text_test[10000:11000]
    order = 5
else:
    print("Slow training-----------------------")
    text_markov = text_markov
    text_test = text_test[10000:20000]
    order = 5

mark = markov.Markov(text_markov, order=order)
sequence = mark.encode(text_test)
text2 = mark.decode(sequence)
assert text2 == text_test

sequence2 = numpy.random.uniform(0, 1.0, 1000)
print(sequence2)
print(mark.decode(sequence2))
