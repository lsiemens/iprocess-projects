import markov
import glob
import numpy
from matplotlib import pyplot

path = "./Books/*.txt"
path2 = "./Books_secondary/*.txt"
replace = [("ë", "e"), ("é", "e"), ("ô","o"), ("à", "a"), ("â", "a"), ("ê", "e"), ("ç", "c"), ("è", "e"), ("ü", "u"), ("ë", "e"), ("æ", "ae"), ("ï", "i")]

text = ""
for file in glob.glob(path):
    booktext = ""
    with open(file, 'r') as fin:
        booktext = fin.read().lower()
        fin.close()
    text += booktext

for (old, new) in replace:
    text = text.replace(old, new)

for char in set(text) - set(markov.Markov.symbols):
    text = text.replace(char, "")


gen_0 = markov.Markov(text, order=0)

for file in glob.glob(path2):
    text = ""
    with open(file, 'r') as fin:
        text = fin.read().lower()
        fin.close()

    for (old, new) in replace:
        text = text.replace(old, new)

    for char in set(text) - set(markov.Markov.symbols):
        text = text.replace(char, "")

    width = 10
    index = numpy.random.random_integers(width, len(text) - width)
    print(index)
    text = text[index - width:index + width]

    print(file)
    print(text)
    sequence_0 = gen_0.getSequence(text)
    print("order 0")
    print(sequence_0, len(sequence_0))

    print(gen_0.getText(sequence_0))
    print("\n\n")
    break
