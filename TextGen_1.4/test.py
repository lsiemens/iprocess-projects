import markov
import huffman
import glob
import numpy
from matplotlib import pyplot

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
sequence = gen_5.getSequence(text_huffman[:10000])
#sequence = gen_5.getSequence(text_huffman)

print(list(set(sequence)))
huf = huffman.Huffman(sequence, list(range(len(gen_5.symbols) + 1)))


sequence2 = gen_5.getSequence(text_huffman[10000:20000])
print("encode")
data = huf.encode(sequence2)
#print(data)
print("decode")
#print(huf.decode(data))
print(sequence2 == huf.decode(data))
