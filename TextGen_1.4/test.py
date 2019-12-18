import markov
import huffman
import utils
import marsta

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

text_markov = "\n".join([texts[file] for file in files[:-5]])
text_huffman = "\n".join([texts[file] for file in files[-5:-3]])
text_huffman_second = "\n".join([texts[file] for file in files[-3:-1]])
text_test = "\n".join([texts[file] for file in [files[-1]]])

if isfast:
    print("Fast training-----------------------")
    text_markov = text_markov[10000:110000]
    text_huffman = text_huffman[10000:20000]
    text_huffman_second = text_huffman_second[10000:20000]
    text_test = text_test[10000:11000]
    width = 4
    order = 5
    num = 100
else:
    print("Slow training-----------------------")
    text_markov = text_markov
    text_huffman = text_huffman
    text_huffman_second = text_huffman_second[10000:500000]
    text_test = text_test[10000:20000]
    width = 8
    order = 5
    num = 10000

mark = markov.Markov(text_markov, order=order)

training_sequence = mark.getSequence(text_huffman)
huff = huffman.Huffman(training_sequence, list(range(len(mark.symbols) + 1)))

training_sequence = mark.getSequence(text_huffman_second)
training_data = huff.encode(training_sequence)
training_data = [training_data[i:i + width] for i in range(0, width*(len(training_data)//width), width)]
huff_second = huffman.Huffman(training_data, [''.join(bit) for bit in itertools.product(*(["01"]*width))])

sequence = mark.getSequence(text_test)
data = huff.encode(sequence)

data = [data[i:i + width] for i in range(0, width*(len(data)//width), width)]
data_second = huff_second.encode(data)
data2 = huff_second.decode(data_second)
data2 = "".join(data2)

sequence2 = huff.decode(data2)
text2 = mark.getText(sequence2)

data = "".join(data)

#print(text)
_txt = ""
for i in range(500):
    if text_test[i] == " ":
        _txt += "||"

    _txt += repr(text_test[i]) + ", " + str(sequence[i]) + "|"

    if text_test[i] == " ":
        _txt += "||"
print(_txt)
#print(sequence)
#print(data)
#print(data_second)
#print("decomp")
#print(data2)
#print(sequence2)
#print(text2)
assert (text2 in text_test)
assert (",".join(map(str,sequence2)) in ",".join(map(str,sequence)))
assert (data2 in "".join(data))

## ---------------------CLEAN
print("test data historical dependence with MarSta")
print("sequence")
tester = marsta.MarSta(sequence, order=1)
tester_d = tester.allIndependent()
for t in tester_d:
    print(t)

print("data")
tester = marsta.MarSta(data, order=3)
tester_d = tester.allIndependent()
for t in tester_d:
    print(t)

print("data_second")
tester = marsta.MarSta(data_second, order=3)
tester_d = tester.allIndependent()
for t in tester_d:
    print(t)
## ---------------------CLEAN

print("data_second ", len(data_second), "bits.")
print("data ", len(data2), "bits.")
print("text ", len(text2)*numpy.log(len(mark.symbols))/numpy.log(2), "bits.")
print("compression data_second ", len(data_second) / (len(text2)*numpy.log(len(mark.symbols))/numpy.log(2)), "percent")
print("compression data ", len(data2) / (len(text2)*numpy.log(len(mark.symbols))/numpy.log(2)), "percent")

tdata = [True if symbol=="1" else False for symbol in data]
print("data moments")

while len(tdata) > 100:
    print(numpy.array(tdata).mean(), numpy.array(tdata).std()/numpy.sqrt(len(tdata) - 1))
    tmp = []
    for i in range(len(tdata)//2):
        tmp.append(tdata[i*2]^tdata[1 + i*2])
    tdata = tmp

tdata = [True if symbol=="1" else False for symbol in data_second]
print("data_second moments")

while len(tdata) > 100:
    print(numpy.array(tdata).mean(), numpy.array(tdata).std()/numpy.sqrt(len(tdata) - 1))
    tmp = []
    for i in range(len(tdata)//2):
        tmp.append(tdata[i*2]^tdata[1 + i*2])
    tdata = tmp

print("\n\n\nTesting")
import sys
data = []
sizes = []
for i in range(num):
    if i % (num//100) == 0:
        sys.stdout.write("\r%d%%" % int(100*i/num))
        sys.stdout.flush()

    size = 500
    new_data_second = "".join(map(str,numpy.random.random_integers(0, 1, size)))
    new_data = huff_second.decode(new_data_second)
    new_data = "".join(new_data)

    new_sequence = huff.decode(new_data)
    new_text = mark.getText(new_sequence)

    data.append((len(new_text), new_text, new_data_second, new_data, new_sequence))
    sizes.append(len(new_text))

sizes.sort(reverse=True)
data.sort(key=lambda node: node[0], reverse=True)
print(numpy.array(sizes))

print(data[0][1], ":text\n------------------------------\n")
print(data[0][4], ":sequence\n")

_txt = ""
for i in range(len(data[0][1])):
    if data[0][1][i] == " ":
        _txt += "||"
    _txt += data[0][1][i] + ", " + str(data[0][4][i]) + "|"
    if data[0][1][i] == " ":
        _txt += "||"
print(_txt)

print("bad\n")

print(data[-1][1], ":text\n------------------------------\n")
print(data[-1][4], ":sequence\n")

_txt = ""
for i in range(len(data[-1][1])):
    if data[-1][1][i] == " ":
        _txt += "||"
    _txt += data[-1][1][i] + ", " + str(data[-1][4][i]) + "|"
    if data[-1][1][i] == " ":
        _txt += "||"
print(_txt)
