#!/usr/bin/env python3

import os
import text_gen

class huffman_compression_static:
    """ 
    Compress a string of integer "symbols" as binary data using huffman coding.
    """
    def compress_array(self, array):
        
    
    def decompress_array(self, data):
        pass

    def _generate_code(self):
        pass

class huffman_compression_dynamic:
    """ 
    Simular to 'huffman_compression_statis' except code is customized for each context in the markov chain
    """
    pass

class markov_encoding(text_gen.markov_chain):
    """ 
    Encode text as a string of integers using a markov chain.
    """
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
                context = context[1:] + text[-1]
            else:
                context = ""
            encoding = encoding[1:]
        return text


    def _encode_char(self, ngram, char):
        index_offset = 0
        ordered_chars = []
        
        while len(ngram) >= 0:
            ordering = []
            for character in self.valid_chars:
                if character not in ordered_chars:
                    if ngram + character in self._transitions:
                        ordering.append((character, self._transitions[ngram + character]))
                    else:
                        ordering.append((character, 0))
            ordering = sorted(ordering, key=lambda value: value[1], reverse=True)
            
            for i, (character, value) in enumerate(ordering):
                if value == 0:
                    index_offset += i
                    break
                else:
                    ordered_chars.append(character)
                
                if character == char:
                    return i + index_offset
            
            if len(ngram) == 0:
                raise ValueError("No matching ngram found.")
            
            ngram = ngram[1:]


    def _decode_char(self, ngram, encoding):
        index_offset = 0
        ordered_chars = []
        
        while len(ngram) >= 0:
            ordering = []
            for character in self.valid_chars:
                if character not in ordered_chars:
                    if ngram + character in self._transitions:
                        ordering.append((character, self._transitions[ngram + character]))
                    else:
                        ordering.append((character, 0))
            ordering = sorted(ordering, key=lambda value: value[1], reverse=True)
            
            for i, (character, value) in enumerate(ordering):
                if value == 0:
                    index_offset += i
                    break
                else:
                    ordered_chars.append(character)
                
                if i + index_offset == encoding:
                    return character
            
            if len(ngram) == 0:
                raise ValueError("No matching ngram found.")
            
            ngram = ngram[1:]


    def _check_compatibility(self):
        if not self.lower_order:
            raise ValueError("Transitions not compattible with encoding, recompute with \'lower_order\' enabled")



text = "hello dad,\nthis is a test message. what are you up to today? i am doing well,\n\nluke"

data = markov_encoding()
data.load("sav.dat")
seed, encoding = data.encode_text(text)
print(encoding)
msg = data.decode_text(seed, encoding)
print(msg)
