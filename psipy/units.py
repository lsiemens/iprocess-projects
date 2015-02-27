"""

Units is a package to automaticaly handle the conversion of units for
quantum mechanics solvers defined in the Hartree system of units.

When using this package hbar = me = e = L = 1, where the SI units value
of L is given when the classes are instansiated. All values used in
equations should be defined as unitless multiples of the constants hbar,
me, e and L.

"""

import numpy

_hbar = 1.05457173*10**(-34)#kg m^2 / s
_me = 9.10938291*10**(-31)#kg
_e = 1.60217657*10**(-19)#coulombs

class _value:
    """
    When instanciated _value will represent one fundemental quantity.
    If the mode is "rel" values are displayed as Hartree units if the
    mode is "abs" the values are displayed as SI units.
    
    """
    def __init__(self, value=1, scale=1, rel="", abs="", mode="rel", format=""):
        self._value = value
        self._scale = scale
        self._rel = rel
        self._abs = abs
        self._mode = mode
        self._format = format
        if not (self._mode in ["rel", "abs"]):
            raise ValueError("The mode " + str(self._mode) + " is not a valid mode.")
    
    def __mul__(self, other):
        try:
            if other._value:
                raise TypeError("two instances of cass \"value\" cannot be multiplied.")
        except AttributeError:
            return _value(other*self._value, self._scale, self._rel, self._abs, self._mode, self._format)

    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __str__(self):
        if self._mode == "rel":
            return self._format.format(self._value) + str(self._rel)
        if self._mode == "abs":
            return self._format.format(self._value*self._scale) + str(self._abs)
    
    def get_format(self):
        return self._format

    def set_format(self, format):
        self._format = format

class units:
    def __init__(self, l=1, mode="rel", format="{:e}"):
        #length in meters
        self.l = l
        self.hbar = _hbar
        self.me = _me
        self.e = _e
        self._format = format
        self._mode = mode
        self._E = _value()
        self._F = _value()
        self._P = _value()
        self._T = _value()
        self._M = _value()
        self._L = _value()

    def set_mode(self, mode="rel"):
        self._mode = mode
        self._E._mode = mode
        self._F._mode = mode
        self._P._mode = mode
        self._T._mode = mode
        self._M._mode = mode
        self._L._mode = mode
        
    def get_mode(self):
        return self._mode

    def set_format(self, format=""):
        self._format = format
        self._E._format = format
        self._F._format = format
        self._P._format = format
        self._T._format = format
        self._M._format = format
        self._L._format = format
        
    def get_format(self):
        return self._format

    def get_E(self):
        self._E = _value(1, (self.hbar/self.l)**2/(self.me), r" E_0", "J", self._mode, self._format)
        return self._E

    def get_F(self):
        self._F = _value(1, self.hbar**2/(self.l**3*self.me), r" F_0", "N", self._mode, self._format)
        return self._F
        
    def get_P(self):
        self._P = _value(1, (self.hbar/self.l), " P_0", r"\frac{kg m}{s}", self._mode, self._format)
        return self._P
        
    def get_T(self):
        self._T = _value(1, (self.me*self.l**2)/(self.hbar), " T_0", r"s", self._mode, self._format)
        return self._T

    def get_M(self):
        self._M = _value(1, self.me, " m_e", r"kg", self._mode, self._format)
        return self._M
        
    def get_L(self):
        self._L = _value(1, self.l, " L_0", r"m", self._mode, self._format)
        return self._L

"""
unit = units(l=10**-12, format="{:+1.1e}")
E = unit.get_E()
P = unit.get_P()
T = unit.get_T()
M = unit.get_M()
L = unit.get_L()

print 2.0003*E
print 2.0003*P
print 2.0003*T
print 2.0003*M
print 2.0003*L

unit.set_mode("abs")
unit.set_format("{:1.5e}")

print 2.0003*E
print 2.0003*P
print 2.0003*T
print 2.0003*M
print 2.0003*L


print "hbar", unit.hbar
print "me", unit.me
print "l", unit.l
"""
