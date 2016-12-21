#!/usr/bin/python3

import os
import text_gen

dir = "./prepared_text/"
#text = [dir + "pg1268.txt"]
text = [dir + file for file in os.listdir(dir) if file.startswith("pg")]
order = 5
lower_order = True

chain = text_gen.markov_chain(order, lower_order, valid_chars="abcdefghijklmnopqrstuvwxyz,. \n?\'\"")
chain.calculate_transitions(text)
chain.save("sav.dat")
chain = None

print("saved data")

#gen = text_gen.markov_chain()
#gen.load("sav.dat")
#print("loaded data")
#gen.generate_string(10000)

