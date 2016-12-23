#!/usr/bin/python3

import encode
import blowfish
from os import urandom
from operator import xor

encrypt = False
decrypt = True
password = b"bassword"

# get nonce for blowfish in counter mode
nonce = 0 #int.from_bytes(urandom(8), "big")
def enc_counter():
    return blowfish.ctr_counter(nonce, f = xor)
bcipher = blowfish.Cipher(password)

data = encode.markov_encoding()
data.load("transition.data")

if encrypt:
    text = ""
    with open("./pg1257.txt") as fin:
        text = fin.read().lower()[:1000]
    
    text = "".join(char for char in text if char in data.valid_chars)
    with open("original.txt", "w") as fout:
        fout.write(text)
    
    encoding = data.encode_text(text)
    encoding = b"".join(bcipher.encrypt_ctr(encoding.rbytes, enc_counter()))
    with open("compressed.data", "wb") as fout:
        fout.write(encoding)
    encoding = None

if decrypt:
    decoding = None
    with open("compressed.data", "rb") as fin:
        decoding = b"".join(bcipher.decrypt_ctr(fin.read(), enc_counter()))
    decoding = encode.rawbytes(bytearray(decoding))
    msg = data.decode_text(decoding)
    
    with open("decompressed2.txt", "w") as fout:
        fout.write(msg)
