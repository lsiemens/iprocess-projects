import markov
import huffman
import utils

import itertools
import numpy
import random
from matplotlib import pyplot

path = "./Books/*.txt"
texts = utils.open_files(path, markov.Markov.symbols)
files = list(texts.keys())
random.shuffle(files)

text_markov = "\n".join([texts[file] for file in files[:-2]])
text_huffman = "\n".join([texts[file] for file in [files[-2]]])
text_test = "\n".join([texts[file] for file in [files[-1]]])

gen_5 = markov.Markov(text_markov[:100000], order=5)
#gen_5 = markov.Markov(text_markov, order=5)
training_sequence = gen_5.getSequence(text_huffman[:10000])
#training_sequence = gen_5.getSequence(text_huffman)

huf = huffman.Huffman(training_sequence, list(range(len(gen_5.symbols) + 1)))

text = str(text_test[:1000])
sequence = gen_5.getSequence(text)
print("markov encode:", sequence)
data = huf.encode(sequence)
print("huffman encode:", data)
sequence2 = huf.decode(data)
#print("huffman decode:", sequence2)
text2 = gen_5.getText(sequence2)
#print("markov decode:", text2)
#print(huf.decode(data))
assert (text == text2)
assert (sequence2 == sequence)

print("data ", len(data), "bits.")
print("text ", len(text)*numpy.log(len(gen_5.symbols))/numpy.log(2), "bits.")
print("compression ", len(data) / (len(text)*numpy.log(len(gen_5.symbols))/numpy.log(2)), "percent")

tdata = [True if symbol=="1" else False for symbol in data]
print("data moments")

while len(tdata) > 100:
    print(numpy.array(tdata).mean(), numpy.array(tdata).std()/numpy.sqrt(len(tdata) - 1))
    tmp = []
    for i in range(len(tdata)//2):
        tmp.append(tdata[i*2]^tdata[1 + i*2])
    tdata = tmp
