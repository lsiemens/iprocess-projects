import os
import text_gen

class rawbytes:
    """ 
    extended bytearray for working with variable length codes
    """
    def __init__(self, byte_array=None):
        if byte_array is None:
            self.rbytes = bytearray([])
            self.start = 0
            self.end = 0
        else:
            self.rbytes = byte_array
            self.start = 0
            self.end = 8*len(byte_array)
            
    def __len__(self):
        return self.end - self.start
    
    def _sanity_check(self):
        if self.start > self.end:
            raise ValueError("Invalid start and end indices.")

        if self.start >= 8:
            raise ValueError("Leading bytes not cleared.")
        
        if self.end > 8*len(self.rbytes):
            raise ValueError("End index out side of bounds.")
        
        if self.end%8 == 0:
            if self.end//8 != len(self.rbytes):
                raise ValueError("Trailing bytes not cleared.")
        else:
            if self.end//8 + 1 != len(self.rbytes):
                raise ValueError("Trailing bytes not cleared.")
    
    def write(self, bools):
        """ 
        Write bool array as bits.
        """
        
        while len(bools) > 0:
            end = self.end%8
            if end == 0:
                self.rbytes.append(0)
            
            mask = 1 << (7 - end)
            if bools[0]:
                self.rbytes[-1] = self.rbytes[-1] | mask
            bools = bools[1:]
            self.end += 1
        
        self._sanity_check()


    def read(self, n):
        """ 
        Read n bits as bool array
        """
        
        if n > self.end - self.start:
            raise ValueError("Cannot read " + str(n) + " bits, index out of bounds.")
        
        data = []
        while n > 0:
            mask = 1 << (7 - self.start%8)
            if self.rbytes[0]&mask != 0:
                data.append(True)
            else:
                data.append(False)
            n -= 1
            self.start += 1

            if (self.start%8 == 0) and (self.start != 0):
                self.rbytes = self.rbytes[1:]
                self.start -= 8
                self.end -= 8
        
        self._sanity_check()
        return data

class tree:
    def __init__(self, character=None, weight=None, children=None):
        self.character = character
        self.weight = weight
        self.children = children
        self.code = None
        
    def generate_codes(self, prefix=[]):
        self.code = prefix
        if self.children is not None:
            self.children[0].generate_codes(self.code + [False])
            self.children[1].generate_codes(self.code + [True])
        
    def __str__(self):
        text = ""
        if self.character is None:
            text = "char: None"
        else:
            text = "char: " + repr(self.character)
        return text + " weight: " + str(self.weight) + " children: " + str(self.children)

class huffman_compression_static:
    """ 
    Compress a string of integer "symbols" as binary data using huffman coding.
    """
    def __init__(self, markov_chain):
        self.markov_chain = markov_chain
        self.codes = None
    
    def compress_value(self, context, symbol, byte_array):
        """ 
        Write symbol to byte_array compressed using huffman code.
        """
        for value, code in self.codes:
            if symbol == value:
                byte_array.write(code)
                break
        return byte_array

    def decompress_value(self, context, byte_array):
        """ 
        Read symbol from byte_array using huffman code.
        """
        bits = []
        while len(byte_array) > 0:
            bits = bits + byte_array.read(1)

            for value, code in self.codes:
                if bits == code:
                    return value
        raise EOFError("End of data stream.")

    def _initalize_compression(self):
        trees = []
        leafs = []
        for character in self.markov_chain.valid_chars:
            if character in self.markov_chain._transitions:
                trees.append(tree(character, self.markov_chain._transitions[character]))
            else:
                raise ValueError("No frequency data for character: " + repr(character))
        leafs = trees
        
        while len(trees) > 1:
            trees = sorted(trees, key=lambda obj: obj.weight, reverse=True)
            newtree = tree(None, trees[-2].weight + trees[-1].weight, trees[-2:])
            trees = trees[:-2] + [newtree]
        
        huffman_tree = trees[0]
        huffman_tree.generate_codes()
        leafs = sorted(leafs, key=lambda obj: obj.weight, reverse=True)
        self.codes = [[i, leaf.code] for i, leaf in enumerate(leafs)]

class huffman_compression_dynamic:
    """ 
    Simular to 'huffman_compression_statis' except code is customized for each context in the markov chain
    """
    def __init__(self, markov_chain):
        self.markov_chain = markov_chain
        self.codes = None
    
    def compress_value(self, context, symbol, byte_array):
        """ 
        Write symbol to byte_array compressed using huffman code.
        """
        self._generate_code(context)

        for value, code in self.codes:
            if symbol == value:
                byte_array.write(code)
                break
        return byte_array

    def decompress_value(self, context, byte_array):
        """ 
        Read symbol from byte_array using huffman code.
        """
        self._generate_code(context)

        bits = []
        while len(byte_array) > 0:
            bits = bits + byte_array.read(1)

            for value, code in self.codes:
                if bits == code:
                    return value
        raise EOFError("End of data stream.")

    def _initalize_compression(self):
        pass
    
    def _generate_code(self, context):
        trees = []
        leafs = []
        processed_chars = []
        scale = 1
        while len(context) >= 0:
            group_total = 0
            for character in self.markov_chain.valid_chars:
                if character not in processed_chars:
                    if context + character in self.markov_chain._transitions:
                        total += self.markov_chain._transitions[context + character]
            
            if total != 0:
                for character in self.markov_chain.valid_chars:
                    if character not in processed_chars:
                        if context + character in self.markov_chain._transitions:
                            trees.append(tree(character, self.markov_chain._transitions[context + character]/float(scale*group_total)))
                            processed_chars.append(character)
                scale = scale*total

            if len(processed_chars) == len(self.markov_chain.valid_chars):
                break

            if (len(context) == 0) and (len(processed_chars) != len(self.markov_chain.valid_chars)):
                raise ValueError("Character frequency data not found.")
            
            context = context[1:]

        leafs = trees
        
        while len(trees) > 1:
            trees = sorted(trees, key=lambda obj: obj.weight, reverse=True)
            newtree = tree(None, trees[-2].weight + trees[-1].weight, trees[-2:])
            trees = trees[:-2] + [newtree]
        
        huffman_tree = trees[0]
        huffman_tree.generate_codes()
        leafs = sorted(leafs, key=lambda obj: obj.weight, reverse=True)
        self.codes = [[i, leaf.code] for i, leaf in enumerate(leafs)]

class markov_encoding(text_gen.markov_chain):
    """ 
    Encode text as a string of integers using a markov chain.
    """
    def __init__(self, order=0, lower_order=False, valid_chars=None):
        super().__init__(order, lower_order, valid_chars)
        #self.compression = huffman_compression_static(self)
        self.compression = huffman_compression_dynamic(self)

    def encode_text(self, text):
        self._check_compatibility()
        encoding = rawbytes()
        context = ""
        text_len = len(text)
        while (len(text) > 0):
            try:
                if len(text)%(text_len//100) == 0:
                    print(100 - len(text)//(text_len//100))
            except ZeroDivisionError:
                pass
            symbol = self._encode_char(context, text[0])
            self.compression.compress_value(context, symbol, encoding)
            if self.order > 0:
                if len(context) >= self.order:
                    context = context[1:]
                context = context + text[0]
            else:
                context = ""
            text = text[1:]
        return encoding

    def decode_text(self, encoding):
        self._check_compatibility()
        text = ""
        context = ""
        encoding_len = len(encoding)
        while (len(encoding) > 0):
            try:
                if len(encoding)%(encoding_len//100) == 0:
                    print(100 - len(encoding)//(encoding_len//100))
            except ZeroDivisionError:
                pass
            
            try:
                symbol = self.compression.decompress_value(context, encoding)
            except EOFError:
                break
            text = text + self._decode_char(context, symbol)
            if self.order > 0:
                if len(context) >= self.order:
                    context = context[1:]
                context = context + text[-1]
            else:
                context = ""
        return text
        
    def load(self, fname):
        super().load(fname)
        self.compression._initalize_compression()

    def _encode_char(self, ngram, char):
        index_offset = 0
        ordered_chars = []
        
        while len(ngram) >= 0:
            ordering = []
            for character in self.valid_chars:
                if character not in ordered_chars:
                    if ngram + character in self._transitions:
                        ordering.append((character, self._transitions[ngram + character]))
                    else:
                        ordering.append((character, 0))
            ordering = sorted(ordering, key=lambda value: value[1], reverse=True)
            
            for i, (character, value) in enumerate(ordering):
                if value == 0:
                    index_offset += i
                    break
                else:
                    ordered_chars.append(character)
                
                if character == char:
                    return i + index_offset
            
            if len(ngram) == 0:
                raise ValueError("No matching ngram found.")
            
            ngram = ngram[1:]


    def _decode_char(self, ngram, encoding):
        index_offset = 0
        ordered_chars = []
        
        while len(ngram) >= 0:
            ordering = []
            for character in self.valid_chars:
                if character not in ordered_chars:
                    if ngram + character in self._transitions:
                        ordering.append((character, self._transitions[ngram + character]))
                    else:
                        ordering.append((character, 0))
            ordering = sorted(ordering, key=lambda value: value[1], reverse=True)
            
            for i, (character, value) in enumerate(ordering):
                if value == 0:
                    index_offset += i
                    break
                else:
                    ordered_chars.append(character)
                
                if i + index_offset == encoding:
                    return character
            
            if len(ngram) == 0:
                raise ValueError("No matching ngram found.")
            
            ngram = ngram[1:]


    def _check_compatibility(self):
        if not self.lower_order:
            raise ValueError("Transitions not compattible with encoding, recompute with \'lower_order\' enabled")
