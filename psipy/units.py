import numpy

_hbar = 1.05457173*10**(-34)#kg m^2 / s
_me = 9.10938291*10**(-31)#kg
_e = 1.60217657*10**(-19)#coulombs

class _value:
    def __init__(self, value=1, scale=1, rel="", abs="", mode="rel"):
        self._value = value
        self._scale = scale
        self._rel = rel
        self._abs = abs
        self._mode = mode
        if not (self._mode in ["rel", "abs"]):
            raise ValueError("The mode " + str(self._mode) + " is not a valid mode.")
    
    def __mul__(self, other):
        try:
            if other._value:
                raise TypeError("two instances of cass \"value\" cannot be multiplied.")
        except AttributeError:
            return _value(other*self._value, self._scale, self._rel, self._abs, self._mode)

    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __str__(self):
        if self._mode == "rel":
            return str(self._value) + str(self._rel)
        if self._mode == "abs":
            return str(self._value*self._scale) + str(self._abs)

class units:
    def __init__(self, l=1, mode="rel"):
        #length in meters
        self.l = l
        self.hbar = _hbar
        self.me = _me
        self.e = _e
        self._mode = mode
        self._E = _value()
        self._P = _value()
        self._T = _value()
        self._M = _value()
        self._L = _value()

    def set_mode(self, mode="rel"):
        self._mode = mode
        self._E._mode = mode
        self._P._mode = mode
        self._T._mode = mode
        self._M._mode = mode
        self._L._mode = mode
        
    def get_mode(self):
        return self._mode

    def get_E(self):
        self._E = _value(1, (numpy.pi*self.hbar/self.l)**2/(2*self.me), r" E_0", "J", self._mode)
        return self._E
        
    def get_P(self):
        self._P = _value(1, (numpy.pi*self.hbar/self.l), " P_0", r"\frac{kg m}{s}", self._mode)
        return self._P
        
    def get_T(self):
        self._T = _value(1, (4*self.me*self.l**2)/(numpy.pi*self.hbar), " T_0", r"s", self._mode)
        return self._T

    def get_M(self):
        self._M = _value(1, self.me, " m_e", r"kg", self._mode)
        return self._M
        
    def get_L(self):
        self._L = _value(1, self.l, " L_0", r"m", self._mode)
        return self._L

unit = units(l=10**-12)
E = unit.get_E()
P = unit.get_P()
T = unit.get_T()
M = unit.get_M()
L = unit.get_L()

print 2*E
print 3*P
print 4*T
print 5*M
print 6*L

unit.set_mode("abs")

print 1*E
print 1*P
print 1*T
print 1*M
print 1*L

print "hbar", unit.hbar
print "me", unit.me
print "l", unit.l
