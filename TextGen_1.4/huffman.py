class Node:
    def __init__(self, weight):
        self.weight = weight
        self.parent = None
        self.isLarge = None
        self.markov_symbol = None

    def setMarkov(self, symbol):
        self.markov_symbol = symbol
        return self

    def getHuffmanSymbol(self):
        if self.parent is None:
            return ""

        if self.isLarge:
            symbol = "1"
        else:
            symbol = "0"

        return self.parent.getHuffmanSymbol() + symbol

class Huffman:
    def __init__(self, training_data, markov_symbols):
        self._training_data = training_data
        self._markov_symbols = markov_symbols
        self._huffman_symbols = None
        self._encode = {} #self._encode[markov_symbol] == huffman_symbol
        self._decode = {} #self._decode[huffman_symbol] == markov_symbol

        self.train()

    def encode(self, sequence):
        data = ""
        for symbol in sequence:
            data += self._encode[symbol]
        return data

    def decode(self, data):
        sequence = []
        symbol = ""
        while len(data) > 0:
            symbol += data[0]
            data = data[1:]
            if symbol in self._decode:
                sequence.append(self._decode[symbol])
                symbol = ""

            assert len(symbol) <= max(map(len, self._huffman_symbols))

#        print("Remainder:", symbol)
        return sequence

    def train(self):
        frequency_data = {symbol:0 for symbol in self._markov_symbols}
        for key in self._training_data:
            if key in frequency_data:
                frequency_data[key] += 1
            else:
                frequency_data[key] = 1
        frequency_data = list(frequency_data.items())
        frequency_data.sort(key=lambda item: item[0])
        frequency_data.sort(key=lambda item: item[1], reverse=True)
        nodes = [Node(weight).setMarkov(symbol) for (symbol, weight) in frequency_data]
        temp_nodes = list(nodes)

        while len(temp_nodes) > 1:
            temp_nodes = self._combine_nodes(temp_nodes)

        for node in nodes:
            huffman_symbol = node.getHuffmanSymbol()
            self._encode[node.markov_symbol] = huffman_symbol
            self._decode[huffman_symbol] = node.markov_symbol
            print(node.weight, node.parent, node.markov_symbol, huffman_symbol) #--------------------------------

        self._huffman_symbols = [self._encode[symbol] for symbol in self._markov_symbols]

    def _combine_nodes(self, nodes):
        node_small = nodes[-1]
        node_large = nodes[-2]
        node_new = Node(node_small.weight + node_large.weight)

        node_small.parent = node_new
        node_large.parent = node_new
        node_small.isLarge = False
        node_large.isLarge = True
        assert node_large.weight >= node_small.weight

        nodes = nodes[:-2] + [node_new]
        nodes.sort(key=lambda node: node.weight, reverse=True)

        return nodes
