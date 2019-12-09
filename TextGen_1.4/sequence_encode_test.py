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


text2 = ""
for file in glob.glob(path2):
    booktext = ""
    with open(file, 'r') as fin:
        booktext = fin.read().lower()
        fin.close()
    text2 += booktext

for (old, new) in replace:
    text2 = text2.replace(old, new)

for char in set(text2) - set(markov.Markov.symbols):
    text2 = text2.replace(char, "")

gen_0 = markov.Markov(text, order=0)
gen_5 = markov.Markov(text, order=10)
gen2_5 = markov.Markov(text2, order=10)

for file in glob.glob(path2):
    text = ""
    with open(file, 'r') as fin:
        text = fin.read().lower()
        fin.close()

    for (old, new) in replace:
        text = text.replace(old, new)

    for char in set(text) - set(markov.Markov.symbols):
        text = text.replace(char, "")

    width = 500
    index = numpy.random.random_integers(width, len(text) - width)
    print(index)
    text = text[index - width:index + width]

    print(file)
    print(text)
    sequence_0 = gen_0.getSequence(text)
    print("order 0")
    print(sequence_0, len(sequence_0))

    sequence_5 = gen_5.getSequence(text)
    print("order 5")
    print(sequence_5, len(sequence_5))

    print(gen_0.getText(sequence_0))
    assert (text == gen_0.getText(sequence_0))
    assert (gen_5.getText(sequence_5) == gen_0.getText(sequence_0))
    print("\n\n")

    """
    pyplot.hist(sequence_0, alpha=0.5, label='Order: 0')
    pyplot.hist(sequence_5, alpha=0.5, label='Order: 5')
    pyplot.legend(loc='upper right')
    pyplot.show()
    """

    print("break ----------------")
    print(gen2_5.getText(sequence_5))
