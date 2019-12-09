import sys
import numpy

class Markov:
#    symbols = "abcdefghijklmnopqrstuvwxyz0123456789 .,?\'\n"
    symbols = "abcdefghijklmnopqrstuvwxyz0123456789 .,?\n"

    def __init__(self, training_text, order=5):
        self.order = order
        self._training_text = training_text
        self._data = {}
        self._ordered_symbols = self.symbols
        self._break_symbol = len(self._ordered_symbols) #give it an index 1 greater than the least common character

        self.train()

    def train(self):
        self._data = {}
        for i in range(len(self._training_text) - self.order):

            if i % ((len(self._training_text) - self.order)//100) == 0:
                sys.stdout.write("\r%d%%" % int(100*i/(len(self._training_text) - self.order)))
                sys.stdout.flush()

            substring = self._training_text[i:i + self.order + 1]
            for j in range(self.order + 1):
                key = substring[j:]
                if key in self._data:
                    self._data[key] += 1
                else:
                    self._data[key] = 1
        print(", Training complete.\n")

        options = self.getOptions()
        self._ordered_symbols = "".join([key for (key, value) in options])

    def getCanonicalText(self, length, seed=""):
        sequence = [0]*(length - len(seed))
        return self.getText(sequence, seed)

    def getRandomText(self, length, seed=""):
#        sequence = numpy.random.random_integers(low=0, high=len(self._ordered_symbols) - 1, size=length)
        sequence = numpy.random.random_integers(low=0, high=1, size=length)
        print(sequence)
        return self.getText(sequence, seed)

    def getText(self, sequence, seed=""):
        text = seed
        if self.order > 0:
            ngram = seed[-self.order:]
        else:
            ngram = ""
        for index in sequence:
            if index == self._break_symbol:
                ngram = ngram[1:]
                continue
            valid_ngram = False
            while not valid_ngram:
                options = self.getOptions(ngram)
                (char, value) = options[index]
                if value == 0:
                    if len(ngram) > 0:
                        ngram = ngram[1:]
                        continue
                    else:
                        raise ValueError("no options.")

                text += char
                ngram += char
                if self.order > 0:
                    ngram = ngram[-self.order:]
                else:
                    ngram = ""
                valid_ngram = True
        return text

    def getSequence(self, text):
        sequence = []
        ngram = ""
        for i in range(len(text)):
            try:
                if i % (len(text)//100) == 0:
                    sys.stdout.write("\r%d%%" % int(100*i/len(text)))
                    sys.stdout.flush()
            except ZeroDivisionError:
                pass

            valid_ngram = False
            while not valid_ngram:
                char = text[i]
                options = self.getOptions(ngram)
                index = list(map(lambda item: item[0], options)).index(char)
                (_, value) = options[index]
                if value == 0:
                    if len(ngram) > 0:
                        sequence.append(self._break_symbol)
                        ngram = ngram[1:]
                        continue
                    else:
                        raise ValueError("no options.")

                sequence.append(index)
                ngram += char
                if self.order > 0:
                    ngram = ngram[-self.order:]
                else:
                    ngram = ""
                valid_ngram = True
        return sequence

    def getOptions(self, ngram=""):
        options = []
        for char in self._ordered_symbols:
            key = ngram + char
            if key in self._data:
                options.append((char, self._data[key]))
            else:
                options.append((char, 0))
        options.sort(key=lambda item: item[1], reverse=True)
        return options
