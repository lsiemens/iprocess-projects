#!/usr/bin/env python3

import os
import text_gen

class markov_encoding(text_gen.markov_chain):
    def encode_text(self, text):
        self._check_compatibility()
        seed = text[:self.order]
        text = text[len(seed):]
        encoding = []
        context = seed

        while (len(text) > 0):
            encoding.append(self._encode_char(context, text[0]))
                
            if self.order > 0:
                context = context[1:] + text[0]
            else:
                context = ""
            text = text[1:]
        return seed, encoding

    def decode_text(self, seed, encoding):
        self._check_compatibility()
        text = seed
        context = seed

        while (len(encoding) > 0):
            text = text + self._decode_char(context, encoding[0])
                
            if self.order > 0:
                context = context[1:] + text[0]
            else:
                context = ""

            encoding = encoding[1:]
        return text

    def _encode_char(self, ngram, char):
        while len(ngram) >= 0:
            ordering = []
            for character in self.valid_chars:
                if ngram + character in self._transitions:
                    ordering.append((character, self._transitions[ngram + character]))
                else:
                    ordering.append((character, 0))
            ordering = sorted(ordering, key=lambda value: value[1], reverse=True)
            
            for i, (character, value) in enumerate(ordering):
                if character == char:
                    return i
                if value == 0:
                    break
    
            if len(ngram) == 0:
                raise ValueError("No matching ngram found.")
    
            ngram = ngram[1:]

    def _decode_char(self, ngram, encoding):
        ordering = []
        for character in self.valid_chars:
            if ngram + character in self._transitions:
                ordering.append((character, self._transitions[ngram + character]))
            else:
                ordering.append((character, 0))
        ordering = sorted(ordering, key=lambda value: value[1], reverse=True)
        
        for i, (character, value) in enumerate(ordering):
            if i == encoding:
                return character
            if value == 0:
                break

        if len(ngram) == 0:
            raise ValueError("No matching ngram found.")

        return self._decode_char(ngram[1:], encoding)

    def _check_compatibility(self):
        if not self.lower_order:
            raise ValueError("Transitions not compattible with encoding, recompute with \'lower_order\' enabled")

text = "hello dad,\nthis is a test message. what are you up to today? i am doing well,\n\nluke"

data = markov_encoding()
data.load("sav.dat")
print(len(text))
seed, encoding = data.encode_text(text)
print(len(seed) + len(encoding))
print(encoding)
msg = data.decode_text(seed, encoding)
print(msg)
print(len(msg))
