import numpy as np
from ctypes import CDLL, Structure, c_int, c_double, POINTER

class Complex(Structure):
    _fields_ = [('real', c_double), ('imag', c_double)]

# Load the shared library into ctypes
lib = CDLL("lib/norm.so")

# Specify the argument types [Complex*, int]
lib.compute_norm.argtypes = [POINTER(Complex), c_int]

# Specify the return type
lib.compute_norm.restype = c_double

# Create a numpy complex array
complex_array = np.array([3+4j, 1+2j, 5+6j], dtype=np.complex128)

# Convert the numpy complex array into a C-style Complex array
complex_struct_array = (Complex * len(complex_array))(*[Complex(c.real, c.imag) for c in complex_array])

# Compute the norm using the C++ function
norm = lib.compute_norm(complex_struct_array, len(complex_array))

print(f"The norm of {complex_array} is {norm}")
