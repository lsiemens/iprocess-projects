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
        pass

class Huffman:
    def __init__(self, training_data, markov_symbols):
        self._training_data = training_data
        self._markov_symbols = markov_symbols
        self._huffman_symbols = None
        self.encode = {} #self._encode[markov_symbol] == huffman_symbol
        self.decode = {} #self._encode[huffman_symbol] == markov_symbol
        self._tree = None

        self.train()

    def train(self):
        frequency_data = {}
        for key in self._training_data:
            if key in frequency_data:
                frequency_data[key] += 1
            else:
                frequency_data[key] = 1
        frequency_data = list(frequency_data.items())
        frequency_data.sort(key=lambda item: item[0])
        frequency_data.sort(key=lambda item: item[1], reverse=True)
        print(frequency_data)
        nodes = [Node(weight).setMarkov(symbol) for (symbol, weight) in frequency_data]
        temp_nodes = list(nodes)

        while len(temp_nodes) > 1:
            temp_nodes = self._combine_nodes(temp_nodes)

        for node in nodes:
            print(node.weight, node.markov_symbol, node.getHuffmanSymbol())

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
