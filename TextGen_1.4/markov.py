import sys
import numpy

class Markov:
    symbols = "abcdefghijklmnopqrstuvwxyz0123456789 .,?\"\n"

    def __init__(self, training_text, order=5):
        self.order = order
        self._training_text = training_text
        self._data = {}

        self.train()

    def train(self):
        self._data = {symbol:1 for symbol in self.symbols} # force every symbol to occure atleast once in the zeroth order chain
        for i in range(len(self._training_text) - self.order):

            try:
                if i % ((len(self._training_text) - self.order)//100) == 0:
                    sys.stdout.write("\rMarkov training {0:d}%".format(int(100*i/(len(self._training_text) - self.order))))
                    sys.stdout.flush()
            except ZeroDivisionError:
                pass

            substring = self._training_text[i:i + self.order + 1]
            for j in range(self.order + 1):
                key = substring[j:]
                if key in self._data:
                    self._data[key] += 1
                else:
                    self._data[key] = 1
        sys.stdout.write("\n")
        sys.stdout.flush()
        print(num_found)

#    def getCanonicalText(self, length, seed=""):
#        sequence = [0]*(length - len(seed))
#        return self.getText(sequence, seed)

#    def getRandomText(self, length, seed=""):
#        sequence = numpy.random.random_integers(low=0, high=len(self._ordered_symbols) - 1, size=length)
#        sequence = numpy.random.random_integers(low=0, high=1, size=length)
#        print(sequence)
#        return self.getText(sequence, seed)

    def decode(self, sequence):
        text = ""
        ngram = ""
        for i, fraction in enumerate(sequence):
            try:
                if i % (len(sequence)//100) == 0:
                    sys.stdout.write("\rMarkov decode {0:d}%".format(int(100*i/len(sequence))))
                    sys.stdout.flush()
            except ZeroDivisionError:
                pass

            char = self.decode_char(fraction, ngram)

            text += char
            ngram += char
            if self.order > 0:
                ngram = ngram[-self.order:]
            else:
                ngram = ""
        sys.stdout.write("\n")
        sys.stdout.flush()
        return text

    def encode(self, text):
        sequence = []
        ngram = ""
        for i in range(len(text)):
            try:
                if i % (len(text)//100) == 0:
                    sys.stdout.write("\rMarkov encode {0:d}%".format(int(100*i/len(text))))
                    sys.stdout.flush()
            except ZeroDivisionError:
                pass

            char = text[i]
            fraction = self.encode_char(char, ngram)

            sequence.append(fraction)
            ngram += char
            if self.order > 0:
                ngram = ngram[-self.order:]
            else:
                ngram = ""
        sys.stdout.write("\n")
        sys.stdout.flush()
        return sequence

    def encode_char(self, char, ngram=""):
        options = self.getOptions(ngram)
        index = [option[0] for option in options].index(char)
        fractions = [option[1] for option in options]
        lower = sum(fractions[:index])
        higher = lower + fractions[index]
        return numpy.random.uniform(lower, higher, 1)[0]

    def decode_char(self, value, ngram=""):
        options = self.getOptions(ngram)
        lower = 0
        for (char, fraction) in options:
            if value < lower + fraction:
                break
            else:
                lower += fraction
        return char

    def getOptions(self, ngram=""):
        options = []
        missing = []
        for char in self.symbols:
            key = ngram + char
            if key in self._data:
                options.append((char, float(self._data[key])))
            else:
                missing.append(char)

        if len(missing) > 0:
            sub_options = self.getOptions(ngram[1:])
            sub_options = [option for option in sub_options if (option[0] in missing)]
            sub_total = sum([option[1] for option in sub_options])
            sub_options = [(char, 1.0*num/sub_total) for (char, num) in sub_options]
            options = options + sub_options

        total = sum([option[1] for option in options])
        options = [(char, 1.0*num/total) for (char, num) in options]

        options.sort(key=lambda item: item[1], reverse=True)
        return options
