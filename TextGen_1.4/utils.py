import glob
import os

_replacements = [("ë", "e"), ("é", "e"), ("ô","o"), ("à", "a"), ("â", "a"), ("ê", "e"), ("ç", "c"), ("è", "e"), ("ü", "u"), ("ë", "e"), ("æ", "ae"), ("ï", "i"), ("“", "\""), ("”", "\""), ("’", "\""), ("‘", "\""), ("\'", "\"")]

def open_files(path, valid_symbols):
    texts={}
    for file_path in glob.glob(path):
        _, file = os.path.split(file_path)
        text = ""
        with open(file_path, "r") as fin:
            text = fin.read()
            fin.close()
        texts[file] = clean_text(text, valid_symbols)

    assert len(list(texts.keys())) > 0
    return texts

def clean_text(text, valid_symbols):
    text = text.lower()
    for (old, new) in _replacements:
        text = text.replace(old, new)

    for char in set(text) - set(valid_symbols):
        text = text.replace(char, "")
    return text
