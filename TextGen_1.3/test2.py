#!/usr/bin/python3

import encode

text = "pg1257.txt"

data = encode.markov_encoding()
data.load("transition.data")
seed, encoding = data.encode_text(text)
msg = data.decode_text(seed, encoding)
print(msg)
