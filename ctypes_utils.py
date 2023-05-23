from ctypes import Structure, c_double

class Complex(Structure):
    _fields_ = [('real', c_double), ('imag', c_double)]
