#!/usr/bin/python3

import encode
import random

obj = encode.markov_encoding()
obj.load("sav.dat")
context_length = 5
context = ""
failures = 0
i = 0
while True:
    i += 1
    context = "".join([random.choice(obj.valid_chars) for _ in range(context_length)])
    char = random.choice(obj.valid_chars)
    print("context: \'" + repr(context) + "\' char: \'" + repr(char) + "\'")
    encoding = obj._encode_char(context, char)
    decoding = obj._decode_char(context, encoding)
    print("context: \'" + repr(context) + "\' decoding: \'" + repr(decoding) + "\'\n")
    if char != decoding:
        failures += 1
    print(str(failures) + "/" + str(i))
    
