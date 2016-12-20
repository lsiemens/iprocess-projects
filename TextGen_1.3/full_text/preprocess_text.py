#!/usr/bin/python3

import os

source_dir = "./"
target_dir = "../prepared_text/"

files = [file for file in os.listdir(source_dir) if file.startswith("pg")]

for file in files:
    text = ""
    with open(source_dir + file, "r") as fin:
        text = fin.read()
        
    text = text.split("\n")
    start = None
    end = None
    
    for i, line in enumerate(text):
        if line.startswith("***"):
            start = i
            break
            
    text = text[start + 1:]
    
    for i, line in enumerate(text):
        if line.startswith("***"):
            end = i
            break
    
    text = "\n".join(text[:end])
        
    if os.path.isfile(target_dir + file):
        os.remove(target_dir + file)    
    
    with open(target_dir + file, "w") as fout:
        fout.write(text)
    
    print(file)
