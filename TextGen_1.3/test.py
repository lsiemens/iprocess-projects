#!/usr/bin/python3

import os
import text_gen

dir = "./prepared_text/"
text = [dir + file for file in os.listdir(dir) if file.startswith("pg")]
order = 1
lower_order = True

chain = text_gen.markov_chain(order, lower_order, valid_chars="abcdefghijklmnopqrstuvwxyz,. \n")
chain.calculate_transitions(text)
chain.save("transition.data")
