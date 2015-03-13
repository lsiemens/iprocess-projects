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

    def initalize(self, fname, format, sets, set_size, integrator, l, g, dt, theta, theta_d, theta_dd):
        self.fname = fname
        self.format = format
        self.sets = sets
        self.set_size = set_size
        self.integrator = integrator
        self.l = l
        self.g = g
        self.dt = dt
        self.data = None
        
        self.n = len(theta)
        if (len(theta_d) != self.n) or (len(theta_dd) != self.n):
            raise ValueError("theta, theta_d and theta_dd must all be of the same size")
        
        with open(self.fname, 'w') as file_out:
            self._write_int(file_out, self.format)
            self._write_int(file_out, self.sets)
            self._write_int(file_out, self.set_size)
            self._write_int(file_out, self.integrator)
            self._write_int(file_out, self.n)
            self._write_value(file_out, self.l)
            self._write_value(file_out, self.g)
            self._write_value(file_out, self.dt)
            self._write_array(file_out, theta)
            self._write_array(file_out, theta_d)
            self._write_array(file_out, theta_dd)

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

    def _write_value(self, file_out, value, new_line=True):
        if new_line:
            file_out.write('{0:.16E}'.format(value).rjust(26)+"\n")
        else:
            file_out.write('{0:.16E}'.format(value).rjust(26))

    def _write_int(self, file_out, value):
        file_out.write(str(int(value)).rjust(12)+"\n")

    def _write_array(self, file_out, array):
        for i in xrange(self.n):
            if i < (self.n-1):
                self._write_value(file_out, array[i], False)
            else:
                self._write_value(file_out, array[i])
