#!/usr/bin/python3

import encode

text = ""
with open("./pg1257.txt") as fin:
    text = fin.read().lower()

data = encode.markov_encoding()
data.load("transition.data")

text = "".join(char for char in text if char in data.valid_chars)
with open("raw.txt", "w") as fout:
    fout.write(text)

encoding = data.encode_text(text)
with open("compressed.data", "wb") as fout:
    fout.write(encoding.rbytes)
encoding = None

decoding = None
with open("compressed.data", "rb") as fin:
    decoding = encode.rawbytes(bytearray(fin.read()))
msg = data.decode_text(decoding)

with open("raw_out2.txt", "w") as fout:
    fout.write(msg)
