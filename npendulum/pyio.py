"support data io for interactions with pyio.f90"

import numpy

class pyio:
    def __init__(self, fname):
        self.fname = fname
        self.format = 0
        self.sets = 0
        self.set_size = 0
        self.integrator = 0
        self.n = 0
        self.l = 0
        self.g = 0
        self.dt = 0
        self.data = None

    def load(self):
        with open(self.fname, 'r') as file_in:
            line = file_in.readline()
            self.format = int(line)
            line = file_in.readline()
            self.sets = int(line)
            line = file_in.readline()
            self.set_size = int(line)
            line = file_in.readline()
            self.integrator = int(line)
            line = file_in.readline()
            self.n = int(line)

            line = file_in.readline()
            self.l = float(line)
            line = file_in.readline()
            self.g = float(line)
            line = file_in.readline()
            self.dt = float(line)
            
            self.data = []
            for line in file_in:
                line_data = [float(item) for item in line.split()]
                self.data.append(line_data)
            self.data = numpy.array(self.data)            
