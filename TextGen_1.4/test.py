import markov
import huffman
import glob
import numpy
from matplotlib import pyplot

import marsta

path = "./Books/*.txt"
path2 = "./Books_secondary/*.txt"
replace = [("ë", "e"), ("é", "e"), ("ô","o"), ("à", "a"), ("â", "a"), ("ê", "e"), ("ç", "c"), ("è", "e"), ("ü", "u"), ("ë", "e"), ("æ", "ae"), ("ï", "i")]

text_markov = ""
for file in glob.glob(path):
    booktext = ""
    with open(file, 'r') as fin:
        booktext = fin.read().lower()
        fin.close()
    text_markov += booktext

for (old, new) in replace:
    text_markov = text_markov.replace(old, new)

for char in set(text_markov) - set(markov.Markov.symbols):
    text_markov = text_markov.replace(char, "")


text_huffman = ""
for file in glob.glob(path2):
    booktext = ""
    with open(file, 'r') as fin:
        booktext = fin.read().lower()
        fin.close()
    text_huffman += booktext

for (old, new) in replace:
    text_huffman = text_huffman.replace(old, new)

for char in set(text_huffman) - set(markov.Markov.symbols):
    text_huffman = text_huffman.replace(char, "")

#gen_5 = markov.Markov(text_markov[:100000], order=5)
gen_5 = markov.Markov(text_markov, order=5)
sequence = gen_5.getSequence(text_huffman[:len(text_huffman)//2])
#sequence = gen_5.getSequence(text_huffman)

print(list(set(sequence)))
huf = huffman.Huffman(sequence, list(range(len(gen_5.symbols) + 1)))

text = text_huffman[len(text_huffman)//2:len(text_huffman)//2 + 30000]
#text = str(input("Enter text (all lower case):"))
sequence2 = gen_5.getSequence(text)
#print("markov encode:", sequence2)
data = huf.encode(sequence2)
#print("huffman encode:", data)
sequence3 = huf.decode(data)
#print("huffman decode:", sequence3)
text2 = gen_5.getText(sequence3)
#print("markov decode:", text2)
#print(huf.decode(data))
print(text == text2)
print(sequence2 == sequence3)

tester = marsta.MarSta(sequence2, order=1)
tester_d = tester.allIndependent()
for t in tester_d:
    print(t)

tester = marsta.MarSta(data, order=3)
tester_d = tester.allIndependent()
for t in tester_d:
    print(t)

print("data ", len(data), "bits.")
print("text ", len(text)*numpy.log(len(gen_5.symbols))/numpy.log(2), "bits.")
print("compression ", len(data) / (len(text)*numpy.log(len(gen_5.symbols))/numpy.log(2)), "percent")

huf2 = huffman.Huffman(list(text_huffman[:10000]), list(gen_5.symbols))
#huf2 = huffman.Huffman(list(text_huffman), list(gen_5.symbols))
print(huf2._encode)

data2 = huf2.encode(text)
#print(data2)
text3 = "".join(huf2.decode(data2))
#print(text3)
print(text == text3)

tester = marsta.MarSta(data2, order=3)
tester_d = tester.allIndependent()
for t in tester_d:
    print(t)

print("data ", len(data2), "bits.")
print("text ", len(text)*numpy.log(len(gen_5.symbols))/numpy.log(2), "bits.")
print("compression ", len(data2) / (len(text)*numpy.log(len(gen_5.symbols))/numpy.log(2)), "percent")
