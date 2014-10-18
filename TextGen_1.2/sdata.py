import cPickle

class sdata:
    def __init__(self, shape=(1,)):
        self.shape = shape
        self._data = []
        self._indice = []
        
    def save(self, fname):
        with open(fname, 'ab') as data_file:
            pickler = cPickle.Pickler(data_file, -1)
            pickler.dump(self.shape)
            pickler.dump(self._indice)
        for _ in xrange(len(self._indice)):
            with open(fname, 'ab') as data_file:
                pickler = cPickle.Pickler(data_file,-1)
                pickler.dump(self._data[0])
                del self._data[0]
    
    def load(self, unpickler):
        self.shape = unpickler.load()
        self._indice = unpickler.load()
        self._data = [0]*len(self._indice)
        for i in xrange(len(self._indice)):
            self._data[i] = unpickler.load()
        
    def _normalize(self, valid_chars, first=True):
        if len(self.shape) != 1:
            if first:
                print "0.0, %"
            for i, sub_array in enumerate(self._data):
                sub_array._normalize(valid_chars, False)
                if first:
                    percentage = 100*(i + 1)/float(len(valid_chars))
                    print ("%.1f"%percentage) + " %"
            if first:
                print "\n"
        else:
            total = sum(self._data)
            if total !=0:
                self._data = [value / float(total) for value in self._data]
   
    def __len__(self):
        return self.shape[0]
        
    def __str__(self):
        string = ", ".join([str(data) for data in self])
        return "[" + string + "]"
    
    def __getitem__(self, key):
        try:
            if len(key) == 1:
                key = key[0]
        except TypeError:
            pass

        try:
            index = key[0]
            key_pass = key[1:]
            if index < 0:
                index = index % len(self)

            if (index >= len(self)) or (len(key) >  len(self.shape)):
                raise IndexError("sdata index out of range")
            
            if index in self._indice:
                id = self._indice.index(index)
                return self._data[id][key_pass]
            else:
                tmp_shape = self.shape
                while (len(key) != 0):
                    tmp_shape = tmp_shape[1:]
                    index = key[0]
                    key = key[1:]
                    if index < 0:
                        index = index % len(self)

                    if index >= len(self):
                        raise IndexError("sdata index out of range")
                if len(tmp_shape) != 0:
                    return sdata(tmp_shape)
                else:
                    return 0
        except TypeError:
            if key < 0:
                key = key % len(self)

            if key >= len(self):
                raise IndexError("sdata index out of range")

            if key in self._indice:
                id = self._indice.index(key)
                return self._data[id]
            else:
                if len(self.shape) != 1:
                    return sdata(self.shape[1:])
                else:
                    return 0
    
    def __setitem__(self, key, value):
        try:
            if len(key) == 1:
                key = key[0]
        except TypeError:
            pass
        
        try:
            index = key[0]
            key_pass = key[1:]
            if index < 0:
                index = index % len(self)

            if (index >= len(self)) or (len(key) >  len(self.shape)):
                raise IndexError("sdata index out of range")
            
            if index in self._indice:
                id = self._indice.index(index)
                self._data[id][key_pass] = value
            else:
                self._indice.append(index)
                self._data.append(sdata(self.shape[1:]))
                self._data[-1][key_pass] = value
        except TypeError:
            if key < 0:
                key = key % len(self)

            if key >= len(self):
                raise IndexError("sdata index out of range")
            
            if key in self._indice:
                id = self._indice.index(key)
                self._data[id] = value
            else:
                self._indice.append(key)
                self._data.append(value)
        
    def __iter__(self):
        return sdata_iter(self)

class sdata_iter:
    def __init__(self, sdata_obj):
        self.sdata_obj = sdata_obj
        self.index = 0
        
    def __iter__(self):
        return self
        
    def next(self):
        try:
            data = self.sdata_obj[self.index]
            self.index += 1
            return data
        except IndexError:
            raise StopIteration()
