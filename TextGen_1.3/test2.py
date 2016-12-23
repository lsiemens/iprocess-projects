#!/usr/bin/python3

import encode

text = ""
with open("./full_text/pg1257.txt") as fin:
    text = fin.read()[:10000]

data = encode.markov_encoding()
data.load("transition.data")

text = "".join(char for char in text if char in data.valid_chars)
with open("raw.txt", "w") as fout:
    fout.write(text)

encoding = data.encode_text(text)
msg = data.decode_text(encoding)
print(msg)
