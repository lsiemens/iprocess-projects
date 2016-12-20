#!/usr/bin/python3

####
#
# Copyright (c) 2016, Luke Siemens
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright 
# notice, this list of conditions and the following disclaimer in the 
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its 
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
####

"""
TextGen 1.3

This is a tool to generate random text using a markov chain.

"""

import numpy
import pickle
import os

class markov_chain:
    """ 
    This object can process raw text to produce markov chains.
    """
    
    def __init__(self, order=0, lower_order=False, valid_chars=None):
        """ 
        Parameters
        ----------
        order : Integer
            Order of n-grams to use in the markov chian.
        lower_order : Bool
            Should lower order n-grams also be used to resolve dead ends.
            The default is False.
        proportional : Bool
            If true the ngrams will be selected proportionaly to the 
            frequency of there occurence, else the most frequent ngram
            will allways be selected. The default is False.
        valid_chars : String
            A list of valid characters for the ngrams. If None the string
            'abcdefghijklmnopqrstuvwxyz. \n'. The default is None.

        """

        if (order < 0) and (not isinstance(order, int)):
            raise ValueError("Order must be a positive integer.")
        
        self.order = order
        self.lower_order = lower_order
        if valid_chars is None:
            self.valid_chars = "abcdefghijklmnopqrstuvwxyz. \n"
        else:
            self.valid_chars = valid_chars
        
        # The transitions are stored using a python dict. Since the python 
        # dict is implemented as a hash table it automaticaly resizes as more
        # ngrams are added and once processed accessing the data has O(1) time
        # complexity.
        self._transitions = {}
        self._text = None
        
    def calculate_transitions(self, text):
        """ 
        Read in raw text and caluclate markov state transitions.
        
        Parameters
        ----------
        text : String or List
            The text to be procesed, the name of a file to be processed or
            a list of files to be processed.
        """
        
        if isinstance(text, str):
            if os.path.isfile(text):
                with open(text, "r") as txt_file:
                    text = txt_file.read()
            else:
                raise ValueError("No file found \"" + fname + "\"")
        else:
            tmp_text = ""
            for fname in text:
                if os.path.isfile(fname):
                    with open(fname, "r") as txt_file:
                        tmp_text += txt_file.read() + "\n"
                else:
                    raise ValueError("No file found \"" + fname + "\"")
            text = tmp_text
            
        self._text = "".join(char for char in text.lower() if char in self.valid_chars)
        
        if len(self._text) < self.order:
            raise ValueError("insufficient text to generate ngrams.")
        
        if not self.lower_order:
            #process all ngrams
            for i in range(len(self._text) - (self.order)):
                ngram = self._text[i:i + self.order + 1]
                self._add_ngram(ngram)
        else:
            #process characters in the first ngram
            for tmp_order in range(self.order):
                ngram = self._text[:tmp_order + 1]
                while len(ngram) > 0:
                    self._add_ngram(ngram)
                    ngram = ngram[1:]
            
            #process the rest of the text
            for i in range(len(self._text) - (self.order)):
                ngram = self._text[i:i + self.order + 1]
                while len(ngram) > 0:
                    self._add_ngram(ngram)
                    ngram = ngram[1:]
        
    def generate_string(self, length, seed=None):
        """ 
        Generate random text from markov chain.
        
        Parameters
        ----------
        length : Integer
            Length of the text to generate.
        seed : string
            The inital seed string to start the chain. If None the first ngram
            in the processed text will be used. The default is None.
        
        """
        if seed is None:
            seed = self._text[:self.order + 1]
            
        print("seed (" + str(len(seed)) + "): \"" + seed + "\"")
        print("order: " + str(self.order))
        if not self.lower_order:
            if len(seed) < self.order:
                raise ValueError("the length of the seed must be equal to the order")
        text = seed
        while (len(text) < length):
            if self.order > 0:
                context = text[-self.order:]
            else:
                context = ""
            text += self._get_character(context)
        print(text)
        
    def save(self, fname):
        """ 
        Save markov chain to disk
        
        Parameters
        ----------
        fname : String
            Name of file.
        """
        
        with open(fname, "wb") as fout:
            pickle.dump(self, fout)
        
    def load(self, fname):
        """ 
        Load markov chain from disk.
        
        Parameters
        ----------
        fname : String
            Name of file.
        """
        
        with open(fname, "rb") as fin:
            self._assignment(pickle.load(fin))

    def _add_ngram(self, ngram):
        """ 
        Add ngram to transitions
        
        Parameters
        ----------
        ngram : string
            ngram to add to the markov chain transitions
        """
        if ngram in self._transitions:
            self._transitions[ngram] += 1
        else:
            self._transitions[ngram] = 1
            
    def _get_character(self, ngram):
        char = None
        if not self.lower_order:
            char = self._choose_character(ngram)
        else:
            while True:
                char = self._choose_character(ngram)
                if (char is not None) or (len(ngram) == 0):
                    break
                ngram = ngram[1:]
        if char is None:
            raise ValueError("Failed to find ngram.")
        return char

    def _choose_character(self, ngram):
        character = None
        occurrences = numpy.zeros(shape=(len(self.valid_chars)))
        for i, char in enumerate(self.valid_chars):
            if ngram + char in self._transitions:
                occurrences[i] = self._transitions[ngram + char]
        if occurrences.sum() >= 1.0:
            value = numpy.random.random_integers(0, int(occurrences.sum() - 1), 1)[0]
            for i, char in enumerate(self.valid_chars):
                if occurrences[:i + 1].sum() > value:
                    character = char
                    break
        return character

    def _assignment(self, other):
        """ 
        Implement the code "self = other"
        
        Parameters
        ----------
        other : markov_chain
            Instance of markov_chain to be copied into self.
        """

        self.order = other.order
        self.lower_order = other.lower_order
        self.valid_chars = other.valid_chars
        self._transitions = other._transitions
        self._text = other._text


dir = "./prepared_text/"
#text = [dir + "pg1268.txt"]
text = [dir + file for file in os.listdir(dir) if file.startswith("pg")]
order = 10
lower_order = False

text_gen = markov_chain(order, lower_order, valid_chars="abcdefghijklmnopqrstuvwxyz,. \n?\'\"")
text_gen.calculate_transitions(text)
text_gen.save("sav.dat")
text_gen = None

print("saved data")

gen = markov_chain()
gen.load("sav.dat")
print("loaded data")
gen.generate_string(10000)

